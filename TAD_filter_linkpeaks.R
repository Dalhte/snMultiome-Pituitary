# ------------------------------------------------------------
# 04_TAD_linkpeaks_pipeline.R 
# ------------------------------------------------------------
suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(S4Vectors)
  library(rtracklayer)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(grid)
  library(Seurat)
  library(Signac)
  library(Matrix)
  library(matrixStats)
})

# =======================
# CONFIG (edit paths only)
# =======================
CFG <- list(
  # Inputs/outputs (set paths to your files)
  gtf_file              = "genes.gtf",
  linkpeaks_csv         = "filtered_linkpeaks_all.csv",             # from your previous LinkPeaks export
  tss_bed_out           = "gene_TSS_unique_allGenes.bed",
  tss_from_peaks_bed    = "gene_TSS_unique_peaks.bed",              # if you prefer using peaks annotation
  tad_boundaries_bed    = "rat_boundaries_filtered.bed",            # BED of boundaries (not domains)
  linkpeaks_with_tss    = "filtered_linkpeaks_all_with_TSS.csv",
  linkpeaks_no_TAD      = "filtered_linkpeaks_all_with_TSS_noTAD.csv",
  heatmap_out_tiff      = "Heatmap_DAR_associated_to_DEG.tiff"
)

# =========================================
# 0) Build one unique TSS per gene from GTF
# =========================================
# - One TSS per gene:
#   strand "+": min(start of transcript TSS)
#   strand "-": max(start of transcript TSS)
#   mixed: median(start), strand="*"
# - Writes BED6 (0-based)
build_unique_tss_from_gtf <- function(gtf_file, bed_out) {
  message("[TSS/gtf] Reading GTF → TxDb …")
  
  # Prefer txdbmaker when available to silence the deprecation warning
  txdb <- if (requireNamespace("txdbmaker", quietly = TRUE)) {
    txdbmaker::makeTxDbFromGFF(gtf_file, format = "gtf")
  } else {
    GenomicFeatures::makeTxDbFromGFF(gtf_file, format = "gtf")
  }
  
  # Get transcript starts; TSS is 'start' on '+' and 'end' on '-'
  tx_gr  <- GenomicFeatures::transcripts(txdb, columns = c("gene_id", "tx_name"))
  tss_tx <- GenomicRanges::resize(tx_gr, width = 1, fix = "start", ignore.strand = FALSE)
  df     <- as.data.frame(tss_tx, optional = TRUE)  # keep plain data.frame
  
  # Robust per-gene TSS selection
  # - if unique strand '+' → min(start)
  # - if unique strand '-' → max(start)
  # - mixed/unknown       → median(start) and strand="*"
  suppressWarnings({
    df$seqnames <- as.character(df$seqnames)
    df$strand   <- as.character(df$strand)
    df$gene_id  <- as.character(df$gene_id)
  })
  
  library(dplyr)
  tss_gene <- df %>%
    group_by(gene_id, seqnames) %>%
    summarise(
      .groups = "drop",
      strand_set = list(unique(strand[strand %in% c("+","-")])),
      strand_out = if (length(unlist(strand_set)) == 1) unlist(strand_set)[1] else "*",
      TSS        = {
        s <- if (length(unlist(strand_set)) == 1) unlist(strand_set)[1] else NA_character_
        if (identical(s, "+")) min(start, na.rm = TRUE) else
          if (identical(s, "-")) max(start, na.rm = TRUE) else
            stats::median(start, na.rm = TRUE)
      }
    )
  
  # Build BED6 (0-based half-open) with a 'name' column = gene_id
  bed_df <- data.frame(
    chr    = tss_gene$seqnames,
    start  = pmax(0L, as.integer(tss_gene$TSS) - 1L),  # 0-based
    end    = as.integer(tss_gene$TSS),
    name   = tss_gene$gene_id,                         # BED column 4
    score  = 0L,
    strand = tss_gene$strand_out,
    stringsAsFactors = FALSE
  )
  
  # Deduplicate on (chr,start,end,name) without relying on dplyr::distinct
  keep <- !duplicated(bed_df[, c("chr","start","end","name")])
  bed_df <- bed_df[keep, , drop = FALSE]
  
  stopifnot(nrow(bed_df) > 0)
  
  # Write as BED (tab-separated, no header), compliant with BED6
  data.table::fwrite(
    bed_df[, c("chr","start","end","name","score","strand")],
    file = bed_out, sep = "\t", col.names = FALSE, quote = FALSE
  )
  message("[TSS/gtf] BED written: ", normalizePath(bed_out))
  invisible(bed_df)
}


# =======================================================
# 1) Build one unique TSS per gene from peaks annotation
# =======================================================
# - Uses obj[["peaks"]] annotation; TSS = median of annotated TSS per gene
# Build one unique TSS per gene from peaks annotation 
build_unique_tss_from_peaks <- function(seurat_object, bed_out) {
  # --- Preconditions
  stopifnot(inherits(seurat_object, "Seurat"))
  stopifnot("peaks" %in% names(seurat_object@assays))
  ann <- seurat_object@assays$peaks@annotation
  if (is.null(ann) || length(ann) == 0L)
    stop("[TSS/peaks] 'peaks@annotation' is empty or missing.")
  
  # --- Extract columns from GRanges without creating row.names
  #     (avoid as.data.frame(ann) which enforces unique row.names)
  ann_df <- data.frame(
    seqnames  = as.character(GenomeInfoDb::seqnames(ann)),
    start     = GenomicRanges::start(ann),
    end       = GenomicRanges::end(ann),
    strand    = as.character(GenomicRanges::strand(ann)),
    gene_id   = as.character(S4Vectors::mcols(ann)$gene_id),
    gene_name = as.character(S4Vectors::mcols(ann)$gene_name),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  # --- Keep only rows with a gene_id
  ann_df <- ann_df[!is.na(ann_df$gene_id) & ann_df$gene_id != "", , drop = FALSE]
  if (nrow(ann_df) == 0L) stop("[TSS/peaks] No rows with non-empty gene_id.")
  
  # --- Compute TSS per row (start on '+', end on '-', start if unknown)
  ann_df$TSS <- ifelse(
    ann_df$strand == "+", ann_df$start,
    ifelse(ann_df$strand == "-", ann_df$end, ann_df$start)
  )
  
  # --- Aggregate to one TSS per gene:
  #     + strand consensus: '+' if only '+', '-' if only '-', '*' otherwise
  #     + TSS: min on '+', max on '-', median otherwise
  suppressPackageStartupMessages(library(dplyr))
  tss_gene <- ann_df %>%
    mutate(
      seqnames = as.character(seqnames),
      strand   = ifelse(strand %in% c("+","-"), strand, "*")
    ) %>%
    group_by(gene_id, gene_name, seqnames) %>%
    summarise(
      .groups   = "drop",
      strand_set = list(unique(strand[strand %in% c("+","-")])),
      strand_out = if (length(unlist(strand_set)) == 1) unlist(strand_set)[1] else "*",
      TSS        = {
        s <- if (length(unlist(strand_set)) == 1) unlist(strand_set)[1] else NA_character_
        if (identical(s, "+")) min(TSS, na.rm = TRUE) else
          if (identical(s, "-")) max(TSS, na.rm = TRUE) else
            stats::median(TSS, na.rm = TRUE)
      }
    )
  
  # --- Build BED6 (0-based)
  bed_df <- data.frame(
    chr    = tss_gene$seqnames,
    start  = pmax(0L, as.integer(tss_gene$TSS) - 1L),
    end    = as.integer(tss_gene$TSS),
    name   = ifelse(is.na(tss_gene$gene_name) | tss_gene$gene_name == "",
                    tss_gene$gene_id, tss_gene$gene_name),
    score  = 0L,
    strand = tss_gene$strand_out,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  
  # --- Deduplicate exact BED rows if any remain
  keep <- !duplicated(bed_df[, c("chr","start","end","name")])
  bed_df <- bed_df[keep, , drop = FALSE]
  if (nrow(bed_df) == 0L) stop("[TSS/peaks] Empty BED after deduplication.")
  
  # --- Write BED without header
  data.table::fwrite(
    bed_df[, c("chr","start","end","name","score","strand")],
    file = bed_out, sep = "\t", col.names = FALSE, quote = FALSE
  )
  message("[TSS/peaks] BED written: ", normalizePath(bed_out))
  invisible(bed_df)
}

# =====================================
# 2) Join LinkPeaks table with TSS BED
# =====================================
# - Keeps one TSS per gene (first occurrence)
# - Writes a CSV adding TSS_chr/TSS_start/TSS_end/TSS columns
# Join LinkPeaks CSV with a TSS BED (robust typing + minimal harmonization)
join_linkpeaks_with_tss <- function(link_csv, tss_bed, out_csv) {
  message("[join] Reading LinkPeaks CSV: ", link_csv)
  df <- read.csv(link_csv, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("gene","Chrx","start","end") %in% names(df))) {
    stop("[join] Required columns missing in LinkPeaks CSV: gene, Chrx, start, end")
  }
  # Ensure expected types for the join key and coordinates
  df$gene  <- as.character(df$gene)
  df$Chrx  <- as.character(df$Chrx)
  df$start <- as.integer(df$start)
  df$end   <- as.integer(df$end)
  
  message("[join] Importing TSS BED: ", tss_bed)
  bed <- rtracklayer::import(tss_bed, format = "BED")
  bed_df <- data.frame(
    TSS_chr  = as.character(GenomeInfoDb::seqnames(bed)),
    TSS_start= as.integer(GenomicRanges::start(bed)),
    TSS_end  = as.integer(GenomicRanges::end(bed)),
    gene     = as.character(S4Vectors::mcols(bed)$name),  # force character
    stringsAsFactors = FALSE
  )
  
  # One TSS per gene (keep first; BED may contain duplicates)
  bed_unique <- bed_df |>
    dplyr::group_by(.data$gene) |>
    dplyr::slice(1) |>
    dplyr::ungroup()
  
  # Optional: harmonize chr prefix if styles differ
  add_chr  <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))
  drop_chr <- function(x) sub("^chr", "", x)
  # If one side uses 'chr' and the other not, normalize to 'chr' on both
  has_chr_df  <- any(grepl("^chr", df$Chrx))
  has_chr_tss <- any(grepl("^chr", bed_unique$TSS_chr))
  if (has_chr_df && !has_chr_tss) bed_unique$TSS_chr <- add_chr(bed_unique$TSS_chr)
  if (!has_chr_df && has_chr_tss) {
    df$Chrx           <- add_chr(df$Chrx)
    has_chr_df <- TRUE
  }
  
  # Join on 'gene' as character
  merged <- dplyr::left_join(
    df,
    bed_unique,
    by = "gene"
  )
  
  # Basic diagnostics
  n_all   <- nrow(df)
  n_join  <- sum(!is.na(merged$TSS_start))
  message(sprintf("[join] Matched genes: %d / %d (%.1f%%)", n_join, n_all, 100*n_join/n_all))
  
  # Convenience numeric TSS column
  merged$TSS     <- merged$TSS_start
  merged$TSS_chr <- merged$TSS_chr
  
  utils::write.csv(merged, out_csv, row.names = FALSE)
  message("[join] CSV written: ", normalizePath(out_csv))
  invisible(merged)
}

# ============================================
# 3) Remove links that cross a TAD boundary
# ============================================
# - For each link, create the segment between TSS and nearest peak edge
# - Drop links where that segment intersects any boundary
filter_links_by_TAD_boundaries <- function(links_csv_in, tad_bed, links_csv_out) {
  message("[TAD] Reading links: ", links_csv_in)
  df <- read.csv(links_csv_in, stringsAsFactors = FALSE)
  stopifnot(all(c("gene","Chrx","start","end","TSS","TSS_chr") %in% colnames(df)))
  
  tad <- fread(tad_bed, header = FALSE)
  stopifnot(ncol(tad) >= 3)
  colnames(tad)[1:3] <- c("chr","start","end")
  tad_gr <- makeGRangesFromDataFrame(tad, seqnames.field = "chr",
                                     start.field = "start", end.field = "end",
                                     ignore.strand = TRUE)
  
  add_chr <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))
  df$Chrx    <- add_chr(as.character(df$Chrx))
  df$TSS_chr <- add_chr(as.character(df$TSS_chr))
  
  df <- df |>
    mutate(
      TSS = as.integer(TSS),
      dist_start = abs(start - TSS),
      dist_end   = abs(end - TSS),
      edge       = ifelse(dist_start < dist_end, start, end)
    )
  
  same_chr <- df$Chrx == df$TSS_chr & !is.na(df$TSS) & !is.na(df$edge)
  
  seg_df <- df[same_chr, c("Chrx","TSS","edge")]
  seg_df$seg_start <- pmin(seg_df$TSS, seg_df$edge)
  seg_df$seg_end   <- pmax(seg_df$TSS, seg_df$edge)
  seg_gr <- GRanges(seqnames = seg_df$Chrx, ranges = IRanges(seg_df$seg_start, seg_df$seg_end))
  
  ov <- findOverlaps(seg_gr, tad_gr, ignore.strand = TRUE)
  crosses <- rep(FALSE, nrow(df))
  if (length(ov) > 0) {
    crosses[which(same_chr)[unique(queryHits(ov))]] <- TRUE
  }
  
  df$cross_TAD_boundary <- crosses
  kept <- dplyr::filter(df, !cross_TAD_boundary)
  
  message(sprintf("[TAD] n=%d, crossing=%d, kept=%d (%.1f%%)",
                  nrow(df), sum(crosses), nrow(kept), 100*nrow(kept)/nrow(df)))
  write.csv(kept, links_csv_out, row.names = FALSE)
  message("[TAD] CSV written: ", normalizePath(links_csv_out))
  invisible(kept)
}

# ========================================================
# 4) Heatmap of DAR associated with DEG (regions × clusters)
# ========================================================
# - Reads the no-TAD CSV, extracts 'peak' regions
# - Pulls accessibility from obj[['peaks']]
# - Averages per cluster, row-wise z-score, plots + TIFF
heatmap_DAR_associated_to_DEG <- function(
    links_csv_in  = CFG$linkpeaks_no_TAD,
    seurat_object = obj,
    assay         = "peaks",
    cluster_col   = "source",      # use meta$source if available; else Idents()
    peak_col      = "peak",
    tiff_out      = CFG$heatmap_out_tiff,
    width         = 9,
    height        = 6,
    res           = 300
) {
  message("[heatmap] Reading links: ", links_csv_in)
  data <- read.csv(links_csv_in, stringsAsFactors = FALSE)
  
  stop_if_empty <- function(x, msg) if (length(x) == 0L) stop(msg, call. = FALSE)
  add_chr  <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))
  dash2col <- function(x) sub("^(?:chr)?([A-Za-z0-9]+)-(\\d+)-(\\d+)$", "chr\\1:\\2-\\3", x)
  col2dash <- function(x) sub("^(?:chr)?([A-Za-z0-9]+):(\\d+)-(\\d+)$", "chr\\1-\\2-\\3", x)
  
  stopifnot(assay %in% names(seurat_object@assays))
  feat_all <- rownames(seurat_object[[assay]])
  stop_if_empty(feat_all, sprintf("No features in assay '%s'.", assay))
  feat_use_colon <- any(grepl(":", feat_all))
  feat_have_chr  <- any(grepl("^chr", feat_all))
  
  if (peak_col %in% names(data)) {
    peaks_in <- as.character(data[[peak_col]])
  } else {
    stop("Column 'peak' not found; provide peak strings or rebuild from Chrx/start/end.")
  }
  
  normalize_to_assay <- function(v) {
    x <- v
    x <- if (feat_use_colon) dash2col(x) else col2dash(x)
    has_chr <- grepl("^chr", x)
    if (feat_have_chr && !all(has_chr)) x <- add_chr(x)
    if (!feat_have_chr && any(has_chr)) x <- sub("^chr", "", x)
    x
  }
  peaks_req <- unique(normalize_to_assay(peaks_in))
  
  slot_use <- if ("data" %in% slotNames(seurat_object[[assay]])) "data" else "counts"
  peaks_av  <- intersect(peaks_req, feat_all)
  peaks_mis <- setdiff(peaks_req,  feat_all)
  stop_if_empty(peaks_av, "None of the requested peaks are present in the assay after harmonization.")
  
  mat <- GetAssayData(seurat_object, assay = assay, slot = slot_use)[peaks_av, , drop = FALSE]
  
  cl_vec <- if (cluster_col %in% colnames(seurat_object@meta.data)) {
    seurat_object@meta.data[[cluster_col]]
  } else {
    Idents(seurat_object)
  }
  cl_vec <- as.factor(cl_vec)
  
  design <- Matrix::sparse.model.matrix(~ 0 + cl_vec)
  colnames(design) <- levels(cl_vec)
  n_per_cl <- Matrix::colSums(design)
  avg_by_cl <- as.matrix(mat %*% design)
  avg_by_cl <- sweep(avg_by_cl, 2, n_per_cl, "/")
  
  rmean <- rowMeans(avg_by_cl, na.rm = TRUE)
  rsd   <- matrixStats::rowSds(avg_by_cl, na.rm = TRUE)
  Z     <- sweep(avg_by_cl, 1, rmean, "-")
  Z     <- sweep(Z, 1, ifelse(rsd == 0, 1, rsd), "/")
  Z[!is.finite(Z)] <- 0
  
  brks <- seq(-4, 4, length.out = 201)
  pal  <- colorRampPalette(c("black", "white", "#1aa14e"))(length(brks) - 1)
  
  hm <- pheatmap(
    mat                      = Z,
    cluster_rows             = TRUE,
    cluster_cols             = TRUE,
    scale                    = "none",
    color                    = pal,
    breaks                   = brks,
    show_rownames            = FALSE,
    fontsize_col             = 9,
    main                     = "Accessibility z-scores (noTAD links): regions × clusters",
    clustering_distance_rows = "correlation",
    clustering_method        = "ward.D2",
    legend                   = TRUE,
    border_color             = NA
  )
  
  tiff(filename = tiff_out, units = "in", width = width, height = height, res = res)
  grid::grid.newpage(); grid::grid.draw(hm$gtable); dev.off()
  message("[heatmap] TIFF written: ", normalizePath(tiff_out))
  
  message(length(peaks_av), " regions used; ", length(peaks_mis), " missing (not plotted).")
  if (length(peaks_mis) > 0) message("Missing examples: ", paste(utils::head(peaks_mis, 10), collapse = ", "))
  invisible(list(mat = Z, peaks_used = peaks_av, plot = hm))
}

# ========================================================
# 5) CoveragePlot utilities per gene (Signac)
# ========================================================
# --- Coverage plots per gene (robust to assay name / default assay) ---
coverage_plots_gene <- function(seurat_obj,
                                gene_id,
                                peaks_assay = NULL,
                                subset_cells = NULL,
                                extend_upstream = 250,
                                extend_downstream = 250,
                                heights_cov = 15,
                                out_dir = NULL,
                                print_each = FALSE) {
  # 1) Pick the chromatin assay
  if (is.null(peaks_assay)) {
    peaks_assay <- if ("peaks" %in% Assays(seurat_obj)) "peaks" else
      if ("ATAC"  %in% Assays(seurat_obj)) "ATAC"  else
        stop("No chromatin assay found (neither 'peaks' nor 'ATAC').")
  }
  stopifnot(peaks_assay %in% Assays(seurat_obj))
  
  # 2) Ensure it's a ChromatinAssay (Signac)
  if (!inherits(seurat_obj[[peaks_assay]], "ChromatinAssay")) {
    stop(sprintf("Assay '%s' is not a ChromatinAssay. Set the correct chromatin assay.", peaks_assay))
  }
  
  # 3) Require existing Links (run LinkPeaks beforehand)
  lk <- tryCatch(Signac::Links(seurat_obj[[peaks_assay]]), error = function(e) NULL)
  if (is.null(lk) || length(lk) == 0) {
    stop("No peak→gene Links found on assay '", peaks_assay,
         "'. Run LinkPeaks first (Block 1) and save/reload the object.")
  }
  
  # 4) Extract peak IDs linked to the requested gene
  g <- as.character(gene_id)
  lk_gene <- lk[S4Vectors::mcols(lk)$gene == g]
  if (length(lk_gene) == 0) stop("No link found for gene '", g, "'.")
  
  peak_ids <- as.character(S4Vectors::mcols(lk_gene)$peak)
  peak_ids <- intersect(peak_ids, rownames(seurat_obj[[peaks_assay]]))
  if (length(peak_ids) == 0) stop("Linked peaks for '", g, "' are not present in assay '", peaks_assay, "'.")
  
  # 5) Optional cell subset
  obj_use <- if (is.null(subset_cells)) seurat_obj else subset(seurat_obj, cells = subset_cells)
  
  # 6) Helper to label UCSC like "chr:start-end" in the plot title
  to_ucsc <- function(pk, up = 250, down = 250) {
    s <- strsplit(pk, "-")[[1]]
    st <- max(1, as.integer(s[2]) - up)
    en <- as.integer(s[3]) + down
    paste0(s[1], ":", st, "-", en)
  }
  
  # 7) Generate CoveragePlot for each linked peak
  plots <- lapply(peak_ids, function(pk) {
    Signac::CoveragePlot(
      object            = obj_use,
      assay             = peaks_assay,
      region            = pk,                # "chr-start-end" is accepted by Signac
      expression.assay  = "SoupXRNA",
      links             = TRUE,
      extend.upstream   = extend_upstream,
      extend.downstream = extend_downstream,
      heights           = heights_cov
    ) + ggplot2::ggtitle(to_ucsc(pk, extend_upstream, extend_downstream))
  })
  names(plots) <- vapply(peak_ids, to_ucsc, FUN.VALUE = character(1),
                         up = extend_upstream, down = extend_downstream)
  
  # 8) Optional export
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    for (nm in names(plots)) {
      file_safe <- gsub("[:]", "_", nm)
      ggplot2::ggsave(
        filename = file.path(out_dir, paste0("coverage_", gene_id, "_", file_safe, ".svg")),
        plot     = plots[[nm]],
        device   = "svg",
        width    = 10, height = 6, units = "in"
      )
    }
    message(length(plots), " SVG files written to ", normalizePath(out_dir))
  }
  
  if (isTRUE(print_each)) for (p in plots) print(p)
  invisible(plots)
}

# =====================
# Minimal run examples
# =====================
# 1) Build TSS (choose one):
 # build_unique_tss_from_gtf(CFG$gtf_file, CFG$tss_bed_out)
 # build_unique_tss_from_peaks(obj, CFG$tss_from_peaks_bed)

# 2) Join links + TSS:
 # join_linkpeaks_with_tss(CFG$linkpeaks_csv, CFG$tss_bed_out, CFG$linkpeaks_with_tss)

# 3) Filter by TAD boundaries:
 # filter_links_by_TAD_boundaries(CFG$linkpeaks_with_tss, CFG$tad_boundaries_bed, CFG$linkpeaks_no_TAD)

#4) Heatmap (DAR associated to DEG):
# heatmap_DAR_associated_to_DEG(
#   links_csv_in  = CFG$linkpeaks_no_TAD,
#   seurat_object = obj,
#   assay         = "peaks",
#   cluster_col   = "source"   # or "cluster_clean" if you prefer your cluster labels
# )

# 5) Coverage plots per gene:
# coverage_plots_gene(obj, gene_id = "Fshb", out_dir = "./coverage_Fshb_svg", print_each = FALSE)
