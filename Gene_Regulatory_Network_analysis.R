# =============================================================================
# Script 12 — snMultiome GRN (clean, runnable, RN7 gene annotation included)
#   - Builds/loads RN7 gene annotation (TxDb + OrgDb) into list(genes=GRanges)
#   - SelectTFs (chromVAR vs RNA, TF order from RNA trajectory)
#   - SelectGenes (Peak→Gene with optional TADs, correlation/FDR filter)
#   - TF–gene correlations along trajectory
#   - GRN assembly (motif ∧ P2G ∧ correlation)
#   - GRN scoring + filtering (FDR + correlation) consistent with source script
#   - Weighted network centralities and Gephi exports
#   - Saves key outputs as .rds
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(data.table)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(ComplexHeatmap)
  library(circlize)
  library(Seurat)
  library(Signac)
  library(igraph)
  library(tidygraph)
  # TxDb / OrgDb for RN7:
  library(TxDb.Rnorvegicus.UCSC.rn7.refGene)
  library(org.Rn.eg.db)
  library(GenomicFeatures)
  library(AnnotationDbi)
  library(grid)
  library(ArchR)
  library(scMEGA)
})

# -----------------------------------------------------------------------------
# 0) Configuration — edit paths if needed
# -----------------------------------------------------------------------------
CFG <- list(
  objG_rds           = "objG.rds",
  geneAnnotation_rds = "geneAnnotationRN7.rds",     # created if missing
  tad_bed            = "TADrn7.bed", # set NULL to skip TAD gating
  trajectory_name    = "Trajectory",
  nbins              = 100L,
  smoothWindow       = 7L,
  out_dir            = "results_grn",
  figs_dir           = "figs_grn"
)
dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(CFG$figs_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1) Utilities
# -----------------------------------------------------------------------------
`%||%` <- function(a, b) if (is.null(a) || (is.atomic(a) && length(a) == 0) || anyNA(a)) b else a
get_mat <- function(x) if (inherits(x, "SummarizedExperiment")) SummarizedExperiment::assay(x) else as.matrix(x)
cap_matrix <- function(m, lim = c(-2,2)) { m[m < lim[1]] <- lim[1]; m[m > lim[2]] <- lim[2]; m }
make_scaled <- function(m) {
  m <- get_mat(m)
  if (nrow(m) == 0) return(m)
  if (nrow(m) == 1) {
    z <- (m - mean(m)) / (sd(as.numeric(m)) %||% 1)
    dim(z) <- c(1, ncol(m)); rownames(z) <- rownames(m); colnames(z) <- colnames(m); z
  } else { Z <- t(scale(t(m))); Z[is.na(Z)] <- 0; Z }
}
safe_GetTrajectory <- function(object, assay, trajectory.name, groupEvery = 1, slot = "data",
                               smoothWindow = 7, log2Norm = FALSE) {
  suppressMessages(
    GetTrajectory(object = object, assay = assay, trajectory.name = trajectory.name,
                  groupEvery = groupEvery, slot = slot, smoothWindow = smoothWindow,
                  log2Norm = log2Norm)
  )
}
safe_TrajectoryHeatmap <- function(traj, varCutOff = 0, maxFeatures = NULL,
                                   pal = paletteContinuous(set = "horizonExtra"),
                                   limits = c(-2,2), name = NULL) {
  suppressMessages(
    TrajectoryHeatmap(traj, varCutOff = varCutOff, maxFeatures = maxFeatures,
                      pal = pal, limits = limits, name = name, returnMatrix = TRUE)
  )
}

get_motif_names <- function(obj, atac.assay = "peaks") {
  mm <- try(obj@assays[[atac.assay]]@motifs@motif.names, silent = TRUE)
  if (inherits(mm, "try-error") || is.null(mm)) mm <- try(Motif(obj, assay = atac.assay)@motif.names, silent = TRUE)
  if (inherits(mm, "try-error") || is.null(mm)) stop("Motif names not found in object.")
  make.unique(stringr::str_to_title(mm))
}

read_tads <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path)) return(NULL)
  if (!file.exists(path)) stop("TAD file not found: ", path)
  
  df <- tryCatch(
    data.table::fread(path, sep = "\t", header = FALSE, data.table = FALSE,
                      quote = "", showProgress = FALSE),
    error = function(e) utils::read.table(path, sep = "\t", header = FALSE,
                                          stringsAsFactors = FALSE, comment.char = "")
  )
  if (nrow(df)) {
    bad <- grepl("^(track|browser|#)", df[[1]])
    if (any(bad)) df <- df[!bad, , drop = FALSE]
  }
  if (ncol(df) < 3) stop("TAD BED must have at least 3 columns: chr, start, end.")
  df <- df[, 1:3, drop = FALSE]
  colnames(df) <- c("chr","start","end")
  df$chr   <- as.character(df$chr)
  df$start <- as.integer(df$start)
  df$end   <- as.integer(df$end)
  df <- df[is.finite(df$start) & is.finite(df$end) & nzchar(df$chr), , drop = FALSE]
  if (!nrow(df)) stop("No valid TAD intervals after filtering.")
  if (!any(grepl("^chr", head(df$chr, 10)))) df$chr <- paste0("chr", df$chr)
  
  GenomicRanges::makeGRangesFromDataFrame(
    df, keep.extra.columns = FALSE,
    seqnames.field = "chr", start.field = "start", end.field = "end"
  )
}

# -----------------------------------------------------------------------------
# 2) RN7 gene annotation (create from TxDb+OrgDb if missing; normalize shape)
# -----------------------------------------------------------------------------
normalize_gene_anno <- function(x) {
  # Accept: path (.rds), GRanges, list(genes=GRanges), ArchR GeneAnnotation-like objects
  obj <- x
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) stop("normalize_gene_anno(): file not found: ", x)
    obj <- readRDS(x)
  }
  
  # If it's the raw GRanges
  if (inherits(obj, "GRanges")) {
    g <- obj
  } else if (is.list(obj)) {
    # ArchR createGeneAnnotation typically returns list with $genes
    nm <- names(obj)
    if (length(nm) && "genes" %in% nm && inherits(obj$genes, "GRanges")) {
      g <- obj$genes
    } else {
      stop("normalize_gene_anno(): list provided but no $genes GRanges.")
    }
  } else {
    stop("normalize_gene_anno(): unsupported object type (expect GRanges or list(genes=GRanges)).")
  }
  
  # Ensure UCSC style and required metadata
  suppressWarnings({
    if (!any(grepl("^chr", as.character(GenomicRanges::seqnames(g))[seq_len(min(10L, length(g)))]))) {
      GenomeInfoDb::seqlevelsStyle(g) <- "UCSC"
    }
  })
  mc <- as.data.frame(mcols(g))
  has <- function(n) n %in% colnames(mc)
  symbol    <- if (has("symbol")) mc$symbol else if (has("SYMBOL")) mc$SYMBOL else if (has("gene_name")) mc$gene_name else if (has("gene_id")) mc$gene_id else paste0("gene_", seq_along(g))
  gene_id   <- if (has("gene_id")) mc$gene_id else if (has("ENTREZID")) mc$ENTREZID else as.character(seq_along(g))
  gene_name <- if (has("gene_name")) mc$gene_name else symbol
  mcols(g)$symbol    <- as.character(symbol)
  mcols(g)$gene_id   <- as.character(gene_id)
  mcols(g)$gene_name <- as.character(gene_name)
  
  list(genes = g)
}

make_gene_annotation_RN7 <- function() {
  # Force UCSC style
  invisible(GenomeInfoDb::seqlevelsStyle(TxDb.Rnorvegicus.UCSC.rn7.refGene))
  txdb <- TxDb.Rnorvegicus.UCSC.rn7.refGene
  g0 <- GenomicFeatures::genes(txdb)  # GRanges keyed by ENTREZID
  eid <- names(g0)
  map <- AnnotationDbi::select(org.Rn.eg.db, keys = eid, keytype = "ENTREZID", columns = "SYMBOL")
  sym <- map$SYMBOL[match(eid, map$ENTREZID)]
  mcols(g0)$gene_id   <- eid
  mcols(g0)$symbol    <- as.character(ifelse(is.na(sym) | sym == "", eid, sym))
  mcols(g0)$gene_name <- mcols(g0)$symbol
  GenomeInfoDb::seqlevelsStyle(g0) <- "UCSC"
  list(genes = g0)
}

load_gene_annotation <- function(rds = CFG$geneAnnotation_rds) {
  # Try to read and validate; if invalid or unreadable, rebuild and overwrite.
  if (file.exists(rds)) {
    ok <- FALSE
    ga <- NULL
    try({
      ga <- normalize_gene_anno(rds)
      ok <- is.list(ga) && "genes" %in% names(ga) && inherits(ga$genes, "GRanges")
    }, silent = TRUE)
    if (ok) {
      message("Loaded gene annotation from RDS: ", rds, " [OK]")
      return(ga)
    } else {
      message("Existing RDS is invalid. Rebuilding RN7 gene annotation…")
    }
  } else {
    message("RDS not found. Building RN7 gene annotation…")
  }
  ga <- make_gene_annotation_RN7()
  saveRDS(ga, rds)
  message("Saved gene annotation to: ", rds)
  ga
}

# -----------------------------------------------------------------------------
# 3) Core loaders
# -----------------------------------------------------------------------------
load_objG <- function(path = CFG$objG_rds) {
  if (!file.exists(path)) stop("objG RDS not found: ", path)
  readRDS(path)
}

# -----------------------------------------------------------------------------
# 4) SelectTFs — enforce TF order by RNA trajectory (logic aligned to source)
# -----------------------------------------------------------------------------
SelectTFs <- function(object,
                      tf.assay = "chromvarPeaks",
                      rna.assay = "SoupXRNA",
                      atac.assay = "peaks",
                      trajectory.name = CFG$trajectory_name,
                      groupEvery = 1,
                      p.cutoff = 0.01,
                      cor.cutoff = 0.3) {
  motif_names <- get_motif_names(object, atac.assay)
  
  trajMM  <- safe_GetTrajectory(object, tf.assay, trajectory.name, groupEvery,
                                slot = "data", smoothWindow = CFG$smoothWindow,
                                log2Norm = FALSE)
  rownames(trajMM) <- motif_names
  
  trajRNA <- safe_GetTrajectory(object, rna.assay, trajectory.name, groupEvery,
                                slot = "data", smoothWindow = CFG$smoothWindow,
                                log2Norm = TRUE)
  
  df.cor <- GetCorrelation(trajMM, trajRNA)
  df.cor <- df.cor[df.cor$adj_p < p.cutoff & df.cor$correlation > cor.cutoff, , drop = FALSE]
  if (nrow(df.cor) == 0) {
    return(list(tfs = df.cor,
                tf_time_points = tibble(gene = character(), time_point = numeric()),
                gene_time_points = tibble(gene = character(), time_point = numeric()),
                heatmap = NULL))
  }
  keep <- df.cor$tfs
  trajMM_sub  <- trajMM[keep, , drop = FALSE]
  trajRNA_sub <- trajRNA[keep, , drop = FALSE]
  
  matRNA_all <- safe_TrajectoryHeatmap(trajRNA, varCutOff = 0, limits = c(-2,2),
                                       pal = paletteContinuous(set = "horizonExtra"),
                                       name = "Gene expression")
  df_gene_tp <- tibble(gene = rownames(matRNA_all),
                       time_point = seq(1, CFG$nbins, length.out = nrow(matRNA_all)))
  rownames(df_gene_tp) <- df_gene_tp$gene
  
  df_tf_tp <- tibble(gene = keep, time_point = NA_real_)
  rownames(df_tf_tp) <- keep
  matched <- intersect(keep, rownames(df_gene_tp))
  if (length(matched)) df_tf_tp[matched, "time_point"] <- df_gene_tp[matched, "time_point"]
  missing <- setdiff(keep, matched)
  if (length(missing)) {
    mat_rna_sub <- get_mat(trajRNA_sub); bins <- seq(1, CFG$nbins, length.out = ncol(mat_rna_sub))
    for (g in missing) {
      if (g %in% rownames(mat_rna_sub)) {
        x <- as.numeric(mat_rna_sub[g, ]); x <- x - min(x, na.rm = TRUE)
        df_tf_tp[g, "time_point"] <- if (sum(x, na.rm = TRUE) > 0) sum(x * bins, na.rm = TRUE) / sum(x, na.rm = TRUE) else mean(bins)
      } else df_tf_tp[g, "time_point"] <- mean(bins)
    }
  }
  df_tf_tp <- df_tf_tp[order(df_tf_tp$time_point), , drop = FALSE]
  ordered_tfs <- df_tf_tp$gene
  
  df.cor <- df.cor[match(ordered_tfs, df.cor$tfs), , drop = FALSE]
  df.cor$time_point <- df_tf_tp$time_point
  
  matTF  <- make_scaled(trajMM_sub)[ordered_tfs, , drop = FALSE] |> cap_matrix()
  matRNA <- make_scaled(trajRNA_sub)[ordered_tfs, , drop = FALSE] |> cap_matrix()
  
  col_tf  <- circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B"))
  col_rna <- circlize::colorRamp2(c(-2, -1, 0, 1, 2), c("#0C0786", "#5B02A3", "#B5367A", "#E56B5D", "#F0F921"))
  
  ht_tf <- Heatmap(matTF, name = "TF activity", col = col_tf, cluster_rows = FALSE, cluster_columns = FALSE,
                   show_row_names = FALSE, show_column_names = FALSE, use_raster = TRUE)
  lab_width <- ComplexHeatmap::max_text_width(rownames(matTF), gp = grid::gpar(fontsize = 7)) + grid::unit(3, "mm")
  lab_anno <- rowAnnotation(TF = anno_text(rownames(matTF), which = "row", gp = grid::gpar(fontsize = 7),
                                           just = "right", location = 0.5, rot = 0), width = lab_width)
  ht_rna <- Heatmap(matRNA, name = "Gene expression", col = col_rna, cluster_rows = FALSE, cluster_columns = FALSE,
                    show_row_names = FALSE, show_column_names = FALSE, use_raster = TRUE)
  ht <- ht_tf + lab_anno + ht_rna
  
  list(tfs = df.cor[df.cor$tfs %in% rownames(matRNA), , drop = FALSE],
       tf_time_points = df_tf_tp[df_tf_tp$gene %in% rownames(matRNA), , drop = FALSE],
       gene_time_points = df_gene_tp,
       heatmap = ht)
}

# -----------------------------------------------------------------------------
# 5) Peak→Gene (robust; requires list(genes=GRanges) with symbol/gene_id)
# -----------------------------------------------------------------------------
PeakToGene <- function(peak.mat, gene.mat, gene_anno, max.dist = 5e7, tad.file = NULL) {
  stopifnot(is.matrix(peak.mat), is.matrix(gene.mat))
  if (!is.list(gene_anno) || !"genes" %in% names(gene_anno))
    stop("gene_anno must be a list with $genes (GRanges).")
  
  genes <- gene_anno$genes
  GenomeInfoDb::seqlevelsStyle(genes) <- "UCSC"
  
  mc  <- S4Vectors::mcols(genes)
  sym <- if ("symbol" %in% colnames(mc)) mc$symbol
  else if ("gene_name" %in% colnames(mc)) mc$gene_name
  else if ("SYMBOL" %in% colnames(mc)) mc$SYMBOL else NULL
  if (is.null(sym)) stop("Gene annotation lacks a symbol/gene_name/SYMBOL column.")
  sym <- as.character(sym)
  
  gene.use <- intersect(sym, rownames(gene.mat))
  if (!length(gene.use)) stop("No overlap between gene symbols in annotation and gene.mat rownames.")
  genes <- genes[match(gene.use, sym)]
  gene.mat <- gene.mat[gene.use, , drop = FALSE]
  
  gene_start <- ifelse(as.character(BiocGenerics::strand(genes)) %in% "+",
                       BiocGenerics::start(genes), BiocGenerics::end(genes))
  genes2 <- GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(genes),
    ranges   = IRanges::IRanges(gene_start, width = 1),
    name     = gene.use,
    gene_id  = if ("gene_id" %in% colnames(mc)) as.character(mc$gene_id[match(gene.use, sym)]) else gene.use,
    strand   = BiocGenerics::strand(genes)
  )
  seRNA <- SummarizedExperiment::SummarizedExperiment(
    assays    = S4Vectors::SimpleList(RNA = gene.mat),
    rowRanges = genes2
  )
  
  p3 <- stringr::str_split_fixed(rownames(peak.mat), "-", 3)
  if (ncol(p3) < 3) stop("peak.mat rownames must be 'chr-start-end'.")
  chr <- p3[,1]; st <- suppressWarnings(as.integer(p3[,2])); en <- suppressWarnings(as.integer(p3[,3]))
  if (!any(grepl("^chr", head(chr, 10)))) chr <- paste0("chr", chr)
  
  peakSet <- GenomicRanges::GRanges(seqnames = chr, ranges = IRanges::IRanges(start = st, end = en))
  seATAC <- SummarizedExperiment::SummarizedExperiment(
    assays    = S4Vectors::SimpleList(ATAC = peak.mat),
    rowRanges = peakSet
  )
  
  rna_resized <- GenomicRanges::resize(SummarizedExperiment::rowRanges(seRNA), width = 2 * max.dist + 1, fix = "center")
  atac_center <- GenomicRanges::resize(SummarizedExperiment::rowRanges(seATAC), width = 1, fix = "center")
  
  ov <- GenomicRanges::findOverlaps(rna_resized, atac_center, ignore.strand = TRUE)
  o  <- as.data.frame(ov)
  if (!nrow(o)) return(data.frame())
  colnames(o) <- c("gene_idx","peak_idx")
  
  rrRNA <- SummarizedExperiment::rowRanges(seRNA)
  rrAT  <- SummarizedExperiment::rowRanges(seATAC)
  
  o$distance <- GenomicRanges::distance(rrRNA[o$gene_idx], rrAT[o$peak_idx])
  o$gene <- S4Vectors::mcols(rrRNA)[o$gene_idx, "name"]
  o$peak <- paste0(as.character(GenomicRanges::seqnames(rrAT))[o$peak_idx], "-",
                   GenomicRanges::start(rrAT)[o$peak_idx], "-",
                   GenomicRanges::end(rrAT)[o$peak_idx])
  
  tads_gr <- NULL
  if (!is.null(tad.file)) {
    tads_gr <- tryCatch(read_tads(tad.file), error = function(e) { warning("TADs ignored: ", conditionMessage(e)); NULL })
    if (!is.null(tads_gr)) {
      GenomeInfoDb::seqlevelsStyle(tads_gr) <- "UCSC"
      gt <- GenomicRanges::findOverlaps(rrRNA, tads_gr)
      pt <- GenomicRanges::findOverlaps(rrAT,  tads_gr)
      gene_to_tad <- data.frame(gene_idx = gt@from, gene_tad = gt@to)
      peak_to_tad <- data.frame(peak_idx = pt@from, peak_tad = pt@to)
      o <- merge(o, gene_to_tad, by = "gene_idx", all.x = TRUE)
      o <- merge(o, peak_to_tad, by = "peak_idx", all.x = TRUE)
      o <- o[o$gene_tad == o$peak_tad & !is.na(o$gene_tad), ]
      if (!nrow(o)) return(data.frame())
    }
  }
  
  AT <- SummarizedExperiment::assay(seATAC)[, , drop = FALSE]
  RN <- SummarizedExperiment::assay(seRNA)[, , drop = FALSE]
  n  <- ncol(AT)
  
  o$correlation <- vapply(seq_len(nrow(o)), function(i)
    suppressWarnings(stats::cor(AT[o$peak_idx[i], ], RN[o$gene_idx[i], ])),
    numeric(1))
  o$t_stat  <- o$correlation / sqrt(pmax(1 - o$correlation^2, .Machine$double.eps) / (n - 2))
  o$p_value <- 2 * stats::pt(-abs(o$t_stat), n - 2)
  o$fdr     <- stats::p.adjust(o$p_value, method = "fdr")
  o[!is.na(o$fdr), ]
}

# -----------------------------------------------------------------------------
# 6) SelectGenes — assemble Peak→Gene on trajectory matrices
# -----------------------------------------------------------------------------
SelectGenes <- function(object,
                        atac.assay = "peaks",
                        rna.assay  = "SoupXRNA",
                        trajectory.name = CFG$trajectory_name,
                        groupEvery = 1,
                        var.cutoff.gene = 0.90,
                        cor.cutoff = 0,
                        fdr.cutoff = 1e-4,
                        distance.cutoff = 0,
                        gene_anno = load_gene_annotation()) {
  trajRNA  <- safe_GetTrajectory(object, rna.assay,  trajectory.name, groupEvery, slot = "data",
                                 smoothWindow = CFG$smoothWindow, log2Norm = TRUE)
  trajATAC <- safe_GetTrajectory(object, atac.assay, trajectory.name, groupEvery, slot = "data",
                                 smoothWindow = CFG$smoothWindow, log2Norm = TRUE)
  
  groupMatRNA  <- safe_TrajectoryHeatmap(trajRNA,  varCutOff = var.cutoff.gene, limits = c(-2,2))
  groupMatATAC <- safe_TrajectoryHeatmap(trajATAC, varCutOff = 0, maxFeatures = nrow(trajATAC),
                                         pal = paletteContinuous(set = "solarExtra"), limits = c(-2,2))
  
  df_atac_tp <- tibble(peak = rownames(groupMatATAC),
                       atac_time_point = seq(1, CFG$nbins, length.out = nrow(groupMatATAC)))
  df_gene_tp <- tibble(gene = rownames(groupMatRNA),
                       gene_time_point = seq(1, CFG$nbins, length.out = nrow(groupMatRNA)))
  
  message("Linking peaks to genes…")
  df.p2g <- PeakToGene(groupMatATAC, groupMatRNA, gene_anno = gene_anno, max.dist = 5e7, tad.file = CFG$tad_bed) %>%
    dplyr::as_tibble() %>%
    dplyr::filter(distance > distance.cutoff, correlation > cor.cutoff, fdr < fdr.cutoff) %>%
    dplyr::left_join(df_atac_tp,  by = "peak") %>%
    dplyr::left_join(df_gene_tp,  by = "gene") %>%
    dplyr::mutate(time_diff = abs(atac_time_point - gene_time_point))
  
  list(p2g = df.p2g, atac_time_points = df_atac_tp, gene_time_points = df_gene_tp)
}

# -----------------------------------------------------------------------------
# 7) TF–Gene correlation and GRN assembly
# -----------------------------------------------------------------------------
GetTFGeneCorrelation <- function(object, tf.use = NULL, gene.use = NULL,
                                 tf.assay = "chromvarPeaks", gene.assay = "SoupXRNA",
                                 atac.assay = "peaks",
                                 trajectory.name = CFG$trajectory_name,
                                 groupEvery = 1) {
  trajMM  <- safe_GetTrajectory(object, tf.assay,  trajectory.name, groupEvery, slot = "data",
                                smoothWindow = CFG$smoothWindow, log2Norm = FALSE)
  trajRNA <- safe_GetTrajectory(object, gene.assay, trajectory.name, groupEvery, slot = "data",
                                smoothWindow = CFG$smoothWindow, log2Norm = TRUE)
  rownames(trajMM) <- get_motif_names(object, atac.assay)
  
  tf_activity     <- safe_TrajectoryHeatmap(trajMM,  varCutOff = 0,    pal = paletteContinuous(set = "solarExtra"))
  gene_expression <- safe_TrajectoryHeatmap(trajRNA, varCutOff = 0.90, pal = paletteContinuous(set = "solarExtra"))
  
  if (!is.null(tf.use))   tf_activity     <- tf_activity[intersect(tf.use,   rownames(tf_activity)), , drop = FALSE]
  if (!is.null(gene.use)) gene_expression <- gene_expression[intersect(gene.use, rownames(gene_expression)), , drop = FALSE]
  
  df.cor <- t(cor(t(tf_activity), t(gene_expression))) %>% as.data.frame()
  if (!is.null(tf.use)) df.cor <- df.cor[, intersect(colnames(df.cor), tf.use), drop = FALSE]
  df.cor$gene <- rownames(df.cor)
  df.cor <- df.cor %>% tidyr::pivot_longer(!gene, names_to = "tf", values_to = "correlation") %>%
    dplyr::select(tf, gene, correlation)
  
  n <- ncol(tf_activity)
  df.cor$t_stat  <- (df.cor$correlation / sqrt(pmax(1 - df.cor$correlation^2, .Machine$double.eps) / (n - 2)))
  df.cor$p_value <- 2 * stats::pt(-abs(df.cor$t_stat), n - 2)
  df.cor$fdr     <- p.adjust(df.cor$p_value, method = "fdr")
  df.cor
}

get_motif_matrix <- function(object, atac.assay = "peaks") {
  mm <- try(object@assays[[atac.assay]]@motifs@data, silent = TRUE)
  if (inherits(mm, "try-error") || is.null(mm)) mm <- try(Motif(object, assay = atac.assay)@data, silent = TRUE)
  if (inherits(mm, "try-error") || is.null(mm)) stop("Motif matrix not found in object.")
  colnames(mm) <- get_motif_names(object, atac.assay)
  mm
}

GetGRN <- function(motif.matching, df.cor, df.p2g, dedup_peaks_tf = TRUE) {
  # Accept base matrix / sparse Matrix / data.frame
  if (is.null(motif.matching)) stop("Provide motif.matching.")
  if (is.data.frame(motif.matching)) motif.matching <- as.matrix(motif.matching)
  if (!inherits(motif.matching, "Matrix") && !is.matrix(motif.matching)) {
    motif.matching <- Matrix::Matrix(as.matrix(motif.matching), sparse = TRUE)
  }
  if (nrow(motif.matching) == 0L || ncol(motif.matching) == 0L) {
    stop("motif.matching is empty after alignment (0 rows or 0 cols).")
  }
  
  # Peak–TF pairs from motif matrix
  summ <- Matrix::summary(Matrix::Matrix(motif.matching * 1, sparse = TRUE))
  df.p2m <- tibble::tibble(
    peak     = rownames(motif.matching)[summ$i],
    tf       = colnames(motif.matching)[summ$j],
    is_bound = as.integer(summ$x)
  )
  if (dedup_peaks_tf) df.p2m <- dplyr::distinct(df.p2m, peak, tf, .keep_all = TRUE)
  
  # Join with peak→gene links (many-to-many is expected)
  df.m2g <- df.p2m |>
    dplyr::inner_join(dplyr::select(df.p2g, peak, gene),
                      by = "peak", relationship = "many-to-many") |>
    dplyr::group_by(tf, gene) |>
    dplyr::summarise(n_peaks = dplyr::n(), .groups = "drop")
  
  # Add TF–gene correlations
  df.grn <- dplyr::left_join(df.m2g, df.cor, by = c("tf", "gene"))
  
  list(df.grn = df.grn, df.m2g = df.m2g)
}

# Simple z-score-based edge scoring (kept for completeness; not used for final GRN anymore)
score_grn <- function(df.grn) {
  stopifnot(all(c("tf","gene","correlation","fdr","n_peaks") %in% names(df.grn)))
  df.grn %>%
    dplyr::mutate(
      regulatory_score = (
        scale(correlation)[,1] +
          scale(log10(n_peaks + 1))[,1] -
          scale(log10(pmax(fdr, 1e-300)))[,1]
      )
    )
}

# -----------------------------------------------------------------------------
# Graph helpers (second, "final" versions, aligned with source script logic)
# -----------------------------------------------------------------------------

# Build directed graph TF -> gene from scored GRN table
build_graph <- function(df.grn_sc) {
  stopifnot(all(c("tf", "gene") %in% names(df.grn_sc)))
  
  edges_all <- df.grn_sc %>%
    dplyr::transmute(
      from = tf,
      to   = gene,
      correlation      = correlation,
      t_stat           = t_stat,
      p_value          = p_value,
      fdr              = fdr,
      n_peaks          = n_peaks,
      regulatory_score = regulatory_score
    )
  
  nodes_all <- tibble::tibble(
    name = unique(c(as.character(edges_all$from),
                    as.character(edges_all$to)))
  )
  
  ig <- igraph::graph_from_data_frame(
    d = edges_all,
    vertices = nodes_all,
    directed = TRUE
  )
  tidygraph::as_tbl_graph(ig)
}

# Add centrality metrics on nodes (weighted, as in source script)
compute_centralities <- function(G) {
  G %>%
    tidygraph::activate(nodes) %>%
    dplyr::mutate(
      degree_in    = tidygraph::centrality_degree(
        mode    = "in",
        weights = abs(regulatory_score)
      ),
      degree_out   = tidygraph::centrality_degree(
        mode    = "out",
        weights = abs(regulatory_score)
      ),
      degree_total = tidygraph::centrality_degree(mode = "total"),
      closeness    = tidygraph::centrality_closeness(),
      betweenness  = tidygraph::centrality_betweenness(
        weights = 1 / pmax(abs(regulatory_score), .Machine$double.eps)
      ),
      eigen        = tidygraph::centrality_eigen(),
      pagerank     = tidygraph::centrality_pagerank(
        weights = abs(regulatory_score)
      )
    )
}

# Attach pseudotime to nodes (from TF and gene time-point tables)
attach_timepoints <- function(G, df_tf_time_point, df_gene_time_point) {
  tf_map <- df_tf_time_point %>%
    dplyr::transmute(
      name     = as.character(gene),
      tf_time  = as.numeric(time_point)
    )
  gene_map <- df_gene_time_point %>%
    dplyr::transmute(
      name      = as.character(gene),
      gene_time = as.numeric(time_point)
    )
  
  time_map <- dplyr::full_join(tf_map, gene_map, by = "name") %>%
    dplyr::mutate(time_point = dplyr::coalesce(tf_time, gene_time)) %>%
    dplyr::select(name, time_point)
  
  G %>%
    tidygraph::activate(nodes) %>%
    dplyr::left_join(time_map, by = "name")
}

# Export nodes/edges for Gephi
export_gephi <- function(G, out_dir, basename = "grn_all") {
  gi <- igraph::as.igraph(G)
  
  nodes_export <- igraph::as_data_frame(gi, what = "vertices") %>%
    dplyr::rename(Id = name) %>%
    dplyr::mutate(
      Id    = as.character(Id),
      Label = Id
    ) %>%
    dplyr::select(
      Id, Label,
      dplyr::any_of(c(
        "time_point",
        "time_phase",
        "degree_in", "degree_out", "degree_total",
        "closeness", "betweenness", "eigen", "pagerank"
      )),
      dplyr::everything()
    ) %>%
    dplyr::distinct(Id, .keep_all = TRUE)
  
  edges_export <- igraph::as_data_frame(gi, what = "edges") %>%
    dplyr::rename(Source = from, Target = to) %>%
    dplyr::mutate(
      Source = as.character(Source),
      Target = as.character(Target)
    )
  
  readr::write_csv(
    nodes_export,
    file.path(out_dir, paste0(basename, "_nodes.csv"))
  )
  readr::write_csv(
    edges_export,
    file.path(out_dir, paste0(basename, "_edges.csv"))
  )
  
  invisible(list(nodes = nodes_export, edges = edges_export))
}

# Helper for gene early/late cut on a 1D score vector
get_cut_idx <- function(x, min_frac = 0.05) {
  n <- length(x)
  if (n < 10) stop("Not enough genes to define early/late split.")
  
  kmin <- max(5L, floor(min_frac * n))
  kmax <- n - kmin
  
  if (requireNamespace("changepoint", quietly = TRUE)) {
    cp <- changepoint::cpt.mean(x, method = "AMOC", penalty = "SIC")
    k  <- changepoint::cpts(cp)[1]
    k  <- min(max(k, kmin), kmax)
  } else {
    ks    <- kmin:kmax
    diffs <- vapply(
      ks,
      function(k) abs(mean(x[1:k]) - mean(x[(k + 1):n])),
      numeric(1)
    )
    k <- ks[which.max(diffs)]
  }
  k
}

# Dot-plot helper for TF metrics
plot_metric <- function(metric,
                        data,
                        log_scale = FALSE,
                        lollipop = FALSE) {
  d <- data %>%
    dplyr::transmute(
      name,
      time_phase,
      value = .data[[metric]],
      tie   = dplyr::coalesce(degree_total, 0)
    ) %>%
    dplyr::filter(is.finite(value))
  
  if (log_scale) {
    d <- d %>%
      dplyr::filter(value > 0) %>%
      dplyr::mutate(value_plot = log10(value))
  } else {
    d <- d %>% dplyr::mutate(value_plot = value)
  }
  
  d <- d %>%
    dplyr::arrange(dplyr::desc(value), dplyr::desc(tie), name)
  d$name <- factor(d$name, levels = rev(d$name))
  
  g <- ggplot2::ggplot(
    d,
    ggplot2::aes(x = value_plot, y = name, color = time_phase)
  )
  if (lollipop) {
    g <- g +
      ggplot2::geom_segment(
        ggplot2::aes(
          x    = min(value_plot, na.rm = TRUE),
          xend = value_plot,
          yend = name
        ),
        linewidth = 0.4
      )
  }
  
  g +
    ggplot2::geom_point(size = 3.5) +
    ggplot2::labs(
      title = paste0(metric, " — TF sorted (high → low)"),
      x     = if (log_scale) paste0(metric, " (log10)") else metric,
      y     = NULL,
      color = "Phase"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "top") +
    ggplot2::scale_color_manual(
      values = c(Early = "#E76F51", Late = "#2A9D8F")
    )
}

# ------------------------------------------------------------------
# 0) Load inputs and ensure output folders
# ------------------------------------------------------------------

objG <- load_objG()
ann  <- load_gene_annotation()  # list(genes = GRanges)

dir.create(CFG$out_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(CFG$figs_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------
# 1) TF selection + Heatmap-1 (TF activity + TF gene expression)
# ------------------------------------------------------------------

selTF <- SelectTFs(objG, cor.cutoff = 0.10)  # p.cutoff = 0.01 by default

if (!is.null(selTF$heatmap)) {
  ComplexHeatmap::draw(
    selTF$heatmap,
    heatmap_legend_side    = "right",
    annotation_legend_side = "right"
  )
  
  svglite::svglite(
    file.path(CFG$figs_dir, "Heatmap1_TFactivity_vs_GeneExpr.svg"),
    width  = 12,
    height = 9
  )
  ComplexHeatmap::draw(
    selTF$heatmap,
    heatmap_legend_side    = "right",
    annotation_legend_side = "right"
  )
  dev.off()
}

df.cor             <- selTF$tfs            # columns: tfs, correlation, adj_p, time_point
df_tf_time_point   <- selTF$tf_time_points # columns: gene, time_point for TFs
df_gene_time_point <- selTF$gene_time_points

saveRDS(df.cor,             file.path(CFG$out_dir, "df_cor_TF_RNA.rds"))
saveRDS(df_tf_time_point,   file.path(CFG$out_dir, "df_tf_time_point.rds"))
saveRDS(df_gene_time_point, file.path(CFG$out_dir, "df_gene_time_point.rds"))

# ------------------------------------------------------------------
# 2) Peak→Gene links on trajectory matrices
# ------------------------------------------------------------------

selGenes <- SelectGenes(
  objG,
  gene_anno       = ann,
  fdr.cutoff      = 1e-4,
  cor.cutoff      = 0,
  distance.cutoff = 0
)
df.p2g <- selGenes$p2g
saveRDS(df.p2g, file.path(CFG$out_dir, "df_p2g.rds"))

# ------------------------------------------------------------------
# 3) TF–Gene correlation on trajectory (restricted to selected TFs + P2G targets)
# ------------------------------------------------------------------

tf.gene.cor <- GetTFGeneCorrelation(
  object        = objG,
  tf.use        = df.cor$tfs,
  gene.use      = unique(df.p2g$gene),
  tf.assay      = "chromvarPeaks",
  gene.assay    = "SoupXRNA",
  trajectory.name = CFG$trajectory_name
)
saveRDS(tf.gene.cor, file.path(CFG$out_dir, "tf_gene_cor.rds"))

# ------------------------------------------------------------------
# 3b) Heatmap-2: TF × target-gene correlation matrix
# ------------------------------------------------------------------

GRNHeatmap <- function(tf.gene.cor,
                       tf.timepoint   = NULL,   # named numeric: TF -> time
                       gene.timepoint = NULL,   # named numeric: gene -> time
                       km = 1) {
  
  stopifnot(all(c("tf","gene","correlation") %in% colnames(tf.gene.cor)))
  
  # 1) Build TF x gene correlation matrix
  mat.cor <- tf.gene.cor %>%
    dplyr::select(tf, gene, correlation) %>%
    tidyr::pivot_wider(names_from = tf, values_from = correlation) %>%
    textshape::column_to_rownames("gene") %>%
    as.matrix()
  
  ## -------------------- COLUMNS = TFs --------------------
  column_ha <- NULL
  if (!is.null(tf.timepoint)) {
    tf_time_sub <- tf.timepoint[colnames(mat.cor)]
    keep_cols   <- !is.na(tf_time_sub)
    mat.cor     <- mat.cor[, keep_cols, drop = FALSE]
    tf_time_sub <- tf_time_sub[keep_cols]
    
    ord_tf      <- order(tf_time_sub)
    mat.cor     <- mat.cor[, ord_tf, drop = FALSE]
    tf_time_sub <- tf_time_sub[ord_tf]
    
    tf_breaks <- seq(min(tf_time_sub, na.rm = TRUE),
                     max(tf_time_sub, na.rm = TRUE),
                     length.out = 100)
    tf_cols   <- ArchR::paletteContinuous(set = "blueYellow", n = 100)
    col_fun   <- circlize::colorRamp2(tf_breaks, tf_cols)
    
    column_ha <- ComplexHeatmap::HeatmapAnnotation(
      tf_time = tf_time_sub,
      col = list(tf_time = col_fun),
      annotation_name_side = "left"
    )
  }
  
  ## -------------------- ROWS = genes ---------------------
  row_ha <- NULL
  if (!is.null(gene.timepoint)) {
    gene_time_sub <- gene.timepoint[rownames(mat.cor)]
    keep_rows     <- !is.na(gene_time_sub)
    mat.cor       <- mat.cor[keep_rows, , drop = FALSE]
    gene_time_sub <- gene_time_sub[keep_rows]
    
    ord_gene      <- order(gene_time_sub)
    mat.cor       <- mat.cor[ord_gene, , drop = FALSE]
    gene_time_sub <- gene_time_sub[ord_gene]
    
    gene_breaks <- seq(min(gene_time_sub, na.rm = TRUE),
                       max(gene_time_sub, na.rm = TRUE),
                       length.out = 100)
    gene_cols   <- ArchR::paletteContinuous(set = "blueYellow", n = 100)
    row_fun     <- circlize::colorRamp2(gene_breaks, gene_cols)
    
    row_ha <- ComplexHeatmap::rowAnnotation(
      gene_time = gene_time_sub,
      col = list(gene_time = row_fun),
      show_legend = TRUE
    )
  }
  
  ## -------------------- Heatmap -------------------------
  ComplexHeatmap::Heatmap(
    mat.cor,
    name = "correlation",
    cluster_columns = FALSE,
    cluster_rows    = is.null(gene.timepoint),
    top_annotation  = column_ha,
    left_annotation = row_ha,
    show_row_names    = FALSE,
    show_column_names = TRUE,
    row_km    = km,
    column_km = km,
    border    = TRUE,
    row_title    = "Genes ordered by pseudotime",
    column_title = "Transcription Factors"
  )
}

if (exists("GRNHeatmap")) {
  gene_time_vec <- setNames(
    df_gene_time_point$time_point,
    df_gene_time_point$gene
  )
  tf_time_vec <- setNames(
    df_tf_time_point$time_point,
    df_tf_time_point$gene
  )
  
  ht2 <- GRNHeatmap(
    tf.gene.cor    = tf.gene.cor,
    tf.timepoint   = tf_time_vec,
    gene.timepoint = gene_time_vec,
    km             = 1
  )
  
  ComplexHeatmap::draw(
    ht2,
    heatmap_legend_side    = "right",
    annotation_legend_side = "right"
  )
  
  svglite::svglite(
    file.path(CFG$figs_dir,
              "Heatmap2_TF_by_TargetGene_correlation.svg"),
    width  = 12,
    height = 10
  )
  ComplexHeatmap::draw(
    ht2,
    heatmap_legend_side    = "right",
    annotation_legend_side = "right"
  )
  dev.off()
}

# ------------------------------------------------------------------
# 4) Motif matrix and GRN assembly
# ------------------------------------------------------------------

mm <- get_motif_matrix(objG, atac.assay = "peaks")
colnames(mm) <- get_motif_names(objG, atac.assay = "peaks")

mm <- mm[
  intersect(rownames(mm), unique(df.p2g$peak)),
  intersect(colnames(mm), unique(tf.gene.cor$tf)),
  drop = FALSE
]

motif.matching <- Matrix::Matrix(mm * 1, sparse = TRUE)
stopifnot(nrow(motif.matching) > 0, ncol(motif.matching) > 0)

cat("motif.matching dims:",
    nrow(motif.matching), "peaks x",
    ncol(motif.matching), "TFs\n")

grn_parts <- GetGRN(motif.matching, tf.gene.cor, df.p2g)
df.grn    <- grn_parts$df.grn
df.m2g    <- grn_parts$df.m2g

saveRDS(df.grn, file.path(CFG$out_dir, "df_grn_raw.rds"))
saveRDS(df.m2g, file.path(CFG$out_dir, "df_m2g_peak2gene.rds"))

# ------------------------------------------------------------------
# 5) GRN filtering + scoring (aligned with source script logic)
# ------------------------------------------------------------------

# 5a) Remove rows with NA numerics
df.grn_filt <- df.grn %>%
  dplyr::filter(dplyr::if_all(dplyr::where(is.numeric), ~ !is.na(.x)))

# 5b) FDR filter on TF–gene correlation (same threshold as source: 1e-4)
df.grn_filtered <- df.grn_filt %>%
  dplyr::filter(fdr < 1e-4)

# 5c) Regulatory score (average of scaled correlation, log10(n_peaks+1), -log10(FDR))
df.grn_sc <- df.grn_filtered %>%
  dplyr::mutate(
    regulatory_score = (
      scale(correlation)[,1] +
        scale(log10(n_peaks + 1))[,1] +
        scale(-log10(fdr + 1e-10))[,1]
    ) / 3
  ) %>%
  # Additional correlation cutoff as in df.grn_final_cor from source script
  dplyr::filter(correlation > 0.75)

saveRDS(df.grn_sc, file.path(CFG$out_dir, "df_grn_scored_filtered.rds"))

# ------------------------------------------------------------------
# 6) Build graph, weighted centralities, pseudotime, and node table
# ------------------------------------------------------------------

G <- build_graph(df.grn_sc)
G <- compute_centralities(G)
G <- attach_timepoints(G, df_tf_time_point, df_gene_time_point)

saveRDS(G, file.path(CFG$out_dir, "GRN_tbl_graph.rds"))

nodes_data_all <- G %>%
  tidygraph::activate(nodes) %>%
  tibble::as_tibble()
readr::write_csv(
  nodes_data_all,
  file.path(CFG$out_dir, "nodes_with_metrics.csv")
)

# ------------------------------------------------------------------
# 7) TF early/late split using pseudotime (Rfx7 / Nfat5 threshold)
# ------------------------------------------------------------------

tf_order <- df_tf_time_point %>%
  dplyr::filter(gene %in% df.cor$tfs)

tp_Rfx7  <- tf_order$time_point[tf_order$gene == "Rfx7"]
tp_Nfat5 <- tf_order$time_point[tf_order$gene == "Nfat5"]

if (length(tp_Rfx7) == 0L || length(tp_Nfat5) == 0L) {
  stop("Rfx7 or Nfat5 not found in df_tf_time_point; cannot define TF early/late split.")
}
cut_tf <- mean(c(tp_Rfx7, tp_Nfat5))

tf_phase_tbl <- tf_order %>%
  dplyr::mutate(
    time_phase = dplyr::if_else(time_point <= cut_tf, "Early", "Late")
  ) %>%
  dplyr::select(gene, time_point, time_phase)

# attach phase to node table
nodes_data_all <- nodes_data_all %>%
  dplyr::left_join(
    tf_phase_tbl %>%
      dplyr::transmute(
        name       = as.character(gene),
        time_phase = factor(time_phase, levels = c("Early", "Late"))
      ),
    by = "name"
  )

readr::write_csv(
  nodes_data_all,
  file.path(CFG$out_dir, "nodes_with_metrics_and_TF_phase.csv")
)

# update G with time_phase
G <- G %>%
  tidygraph::activate(nodes) %>%
  dplyr::left_join(
    nodes_data_all %>% dplyr::select(name, time_phase),
    by = "name"
  )

# ------------------------------------------------------------------
# 8) Gene early/late split (SVD + changepoint on TF–gene correlation)
# ------------------------------------------------------------------

genes_subset <- intersect(
  unique(as.character(df.grn_sc$gene)),
  unique(as.character(tf.gene.cor$gene))
)

mat_gene_tf <- tf.gene.cor %>%
  dplyr::semi_join(tibble::tibble(gene = genes_subset), by = "gene") %>%
  dplyr::select(tf, gene, correlation) %>%
  tidyr::pivot_wider(names_from = tf, values_from = correlation) %>%
  textshape::column_to_rownames("gene") %>%
  as.matrix()

gene_time_map <- setNames(
  as.numeric(df_gene_time_point$time_point),
  as.character(df_gene_time_point$gene)
)
gene_time <- gene_time_map[rownames(mat_gene_tf)]
keep      <- !is.na(gene_time)

mat_gene_tf <- mat_gene_tf[keep, , drop = FALSE]
gene_time   <- gene_time[keep]

ord <- order(gene_time)
mat_gene_tf <- mat_gene_tf[ord, , drop = FALSE]
gene_time   <- gene_time[ord]

mat0 <- scale(mat_gene_tf, center = TRUE, scale = FALSE)
mat0[!is.finite(mat0)] <- 0

sv <- svd(mat0, nu = 1, nv = 0)
gene_score <- as.numeric(sv$u[, 1] * sv$d[1])

k_star  <- get_cut_idx(gene_score)
tp_star <- as.numeric(gene_time[k_star])

cat(sprintf(
  "Gene early/late split: time_point ≈ %.2f | position %d/%d (%.1f%%)\n",
  tp_star, k_star, length(gene_score),
  100 * k_star / length(gene_score)
))

genes_early <- rownames(mat_gene_tf)[seq_len(k_star)]
genes_late  <- rownames(mat_gene_tf)[(k_star + 1):nrow(mat_gene_tf)]

cat(sprintf(
  "Early genes: %d | Late genes: %d\n",
  length(genes_early), length(genes_late)
))

gene_phase_tbl <- tibble::tibble(
  gene       = rownames(mat_gene_tf),
  gene_time  = as.numeric(gene_time),
  gene_phase = dplyr::if_else(
    seq_along(gene_time) <= k_star,
    "Early",
    "Late"
  )
)

saveRDS(
  list(
    tp_star     = tp_star,
    genes_early = genes_early,
    genes_late  = genes_late,
    gene_time   = gene_time,
    gene_phase  = gene_phase_tbl
  ),
  file.path(CFG$out_dir, "gene_early_late_split.rds")
)

# density plot along pseudotime
h <- hist(
  gene_time,
  breaks = 20,
  col    = "grey90",
  border = "grey40",
  main   = "Gene density along pseudotime (GRN subset)",
  xlab   = "pseudotime (1..100)"
)
abline(v = tp_star,           col = "red",    lwd = 2)
abline(v = median(gene_time), col = "grey30", lty = 2)

svg(file.path(CFG$figs_dir, "gene_density_pseudotime.svg"),
    width = 7, height = 5)
plot(
  h,
  col    = "grey90",
  border = "grey40",
  main   = "Gene density along pseudotime (GRN subset)",
  xlab   = "pseudotime (1..100)"
)
abline(v = tp_star,           col = "red",    lwd = 2)
abline(v = median(gene_time), col = "grey30", lty = 2)
dev.off()

# ------------------------------------------------------------------
# 9) Gephi exports + TF network metric dot-plots
# ------------------------------------------------------------------

export_gephi(G, CFG$out_dir, basename = "grn_all")

# TF table with metrics and phase (keep only TFs that were in df.cor$tfs)
df_tf <- nodes_data_all %>%
  dplyr::filter(name %in% unique(as.character(df.cor$tfs))) %>%
  dplyr::mutate(
    time_phase = factor(time_phase, levels = c("Early", "Late"))
  )

metrics_log <- c(
  "degree_in", "degree_out", "degree_total",
  "betweenness", "pagerank", "closeness"
)

p_degree_total <- plot_metric(
  "degree_total",
  data      = df_tf,
  log_scale = TRUE
)
p_degree_total
p_betweenness <- plot_metric(
  "betweenness",
  data      = df_tf,
  log_scale = TRUE
)
p_betweenness
p_closeness <- plot_metric(
  "closeness",
  data      = df_tf,
  log_scale = TRUE
)
p_eigen <- plot_metric(
  "eigen",
  data      = df_tf,
  log_scale = FALSE
)
p_eigen
p_pagerank <- plot_metric(
  "pagerank",
  data      = df_tf,
  log_scale = TRUE
)
p_pagerank
p_closeness <- plot_metric(
  "closeness",
  data      = df_tf,
  log_scale = TRUE
)
p_closeness


ggplot2::ggsave(
  file.path(CFG$figs_dir, "degree_total_TF.svg"),
  p_degree_total, width = 7, height = 9
)
ggplot2::ggsave(
  file.path(CFG$figs_dir, "betweenness_TF.svg"),
  p_betweenness, width = 7, height = 9
)
ggplot2::ggsave(
  file.path(CFG$figs_dir, "closeness_TF.svg"),
  p_closeness, width = 7, height = 9
)
ggplot2::ggsave(
  file.path(CFG$figs_dir, "eigen_TF.svg"),
  p_eigen, width = 7, height = 9
)
ggplot2::ggsave(
  file.path(CFG$figs_dir, "pagerank_TF.svg"),
  p_pagerank, width = 7, height = 9
)
