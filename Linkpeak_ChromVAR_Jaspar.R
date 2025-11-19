# =============================================================================
# Multiome post-processing: Links, Motifs/chromVAR, DAR/DE, checks
# Idempotent, sans accès direct aux slots, exports optionnels.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(ggplot2)
})

set.seed(42)

# ------------------------------ CONFIG ---------------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

CONFIG <- list(
  # Input
  integrated_rds = "C:/Users/Charles/Desktop/Charles/SnMultiome/Data/PGintegrated73.RN7.links.rds",
  
  # Assays
  rna_assay     = "SoupXRNA",
  peaks_pref    = "peaks",
  atac_fallback = "ATAC",
  
  # Genome
  bsgenome_package = "BSgenome.Rnorvegicus.NCBI.rn7.2",
  seqstyle         = "UCSC",
  
  # LinkPeaks
  linkpeaks = list(
    expression_assay = "SoupXRNA",
    distance         = 5e5,
    min_cells        = 10L,
    n_sample         = 200L,
    pvalue_cutoff    = 0.05,
    score_cutoff     = 0.05,
    method           = "pearson",
    genes_use        = NULL
  ),
  
  # chromVAR / motifs
  run_chromvar      = TRUE,
  jaspar_collection = "CORE",
  jaspar_tax_group  = "vertebrates",
  
  # Differential testing
  do_da_by_cluster  = TRUE,
  da_test           = "LR",
  da_latent_var     = NULL,
  da_min_pct        = 0.05,  # DAR: 5%
  
  do_de_by_cluster  = TRUE,  # DEG vs reste
  de_min_pct        = 0,     # DEG: 0%
  
  chromvar_min_pct  = 0.05,  # chromVAR: 5%
  
  # Exports
  do_export = FALSE,
  out_dir   = "data_processed"
)

# ------------------------------ HELPERS --------------------------------------

load_bsgenome <- function(pkg, seqstyle = "UCSC") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("BSgenome '%s' absent. Étapes dépendantes du génome désactivées.", pkg))
    return(NULL)
  }
  g <- get(pkg, envir = asNamespace(pkg))
  GenomeInfoDb::seqlevelsStyle(g) <- seqstyle
  g
}

ensure_groupings <- function(obj) {
  md <- obj@meta.data
  if (!"dataset" %in% colnames(md)) obj$dataset <- obj$orig.ident
  # Respecte ta convention 'cell_Type' si présente
  if ("cell_Type" %in% colnames(md)) {
    obj$cell_Type <- as.character(md$cell_Type)
  } else if (!"cell_type" %in% colnames(md)) {
    obj$cell_type <- as.character(Idents(obj))
  }
  # Identifiant combiné robuste (prend cell_Type si dispo)
  ct_col <- if ("cell_Type" %in% colnames(obj@meta.data)) "cell_Type" else "cell_type"
  obj$grouping_id <- paste0(obj$dataset, "_", obj[[ct_col]])
  obj
}

choose_peaks_assay <- function(obj) {
  assays <- Assays(obj)
  if (CONFIG$peaks_pref %in% assays) return(CONFIG$peaks_pref)
  if (CONFIG$atac_fallback %in% assays) return(CONFIG$atac_fallback)
  stop("Aucun assay chromatin trouvé ('peaks' ou 'ATAC').")
}

run_linkpeaks_once <- function(obj, genome, peaks_assay, cfg) {
  DefaultAssay(obj) <- peaks_assay
  
  existing_links <- tryCatch({
    lk <- Links(obj[[peaks_assay]])
    if (is.null(lk)) 0L else nrow(as.data.frame(lk))
  }, error = function(e) 0L)
  
  if (existing_links > 0) {
    message(sprintf("LinkPeaks: %d liens existants détectés → skip.", existing_links))
    return(obj)
  }
  
  if (!is.null(genome)) {
    main_chroms <- standardChromosomes(genome)
    gr <- granges(obj[[peaks_assay]])
    keep <- as.character(seqnames(gr)) %in% main_chroms
    obj[[peaks_assay]] <- subset(obj[[peaks_assay]], features = rownames(obj[[peaks_assay]])[keep])
    message(sprintf("Filtré: %d pics sur chromosomes standard.", length(which(keep))))
  } else {
    message("Genome NULL → pas de filtre chromosomes / RegionStats.")
  }
  
  if (!is.null(genome)) {
    obj <- RegionStats(obj, genome = genome, assay = peaks_assay)
  }
  
  obj <- LinkPeaks(
    object            = obj,
    peak.assay        = peaks_assay,
    expression.assay  = cfg$expression_assay,
    genes.use         = cfg$genes_use,
    distance          = cfg$distance,
    min.cells         = cfg$min_cells,
    n_sample          = cfg$n_sample,
    pvalue_cutoff     = cfg$pvalue_cutoff,
    score_cutoff      = cfg$score_cutoff,
    method            = cfg$method,
    verbose           = TRUE
  )
  message("LinkPeaks: OK.")
  obj
}

run_chromvar_once <- function(obj, genome, peaks_assay) {
  if (!CONFIG$run_chromvar) { message("chromVAR=FALSE → skip."); return(obj) }
  if (is.null(genome))      { message("Genome NULL → skip motifs/chromVAR."); return(obj) }
  
  has_tfbstools  <- requireNamespace("TFBSTools",  quietly = TRUE)
  has_jaspar2022 <- requireNamespace("JASPAR2022", quietly = TRUE)
  if (!has_tfbstools || !has_jaspar2022) {
    message("TFBSTools/JASPAR2022 manquants → skip motifs/chromVAR.")
    return(obj)
  }
  
  DefaultAssay(obj) <- peaks_assay
  
  jaspar_db <- tryCatch(get("JASPAR2022", envir = asNamespace("JASPAR2022")), error = function(e) NULL)
  if (is.null(jaspar_db)) {
    message("JASPAR2022 non récupérable → skip motifs/chromVAR.")
    return(obj)
  }
  
  if (is.null(Motifs(obj[[peaks_assay]]))) {
    pfm <- TFBSTools::getMatrixSet(
      x    = jaspar_db,
      opts = list(collection = CONFIG$jaspar_collection, tax_group = CONFIG$jaspar_tax_group, all_versions = FALSE)
    )
    obj <- AddMotifs(object = obj, genome = genome, pfm = pfm, assay = peaks_assay)
    message("Motifs ajoutés.")
  } else {
    message("Motifs déjà présents → skip AddMotifs.")
  }
  
  if (!"chromvar" %in% Assays(obj)) {
    obj <- RunChromVAR(object = obj, genome = genome, assay = peaks_assay)
    message("chromVAR: OK (assay 'chromvar').")
  } else {
    message("Assay 'chromvar' déjà présent → skip RunChromVAR.")
  }
  
  obj
}

chromvar_variability_df <- function(obj, assay = "chromvar", layer = "data") {
  if (!requireNamespace("matrixStats", quietly = TRUE)) stop("matrixStats requis.")
  if (!(assay %in% Assays(obj))) stop(sprintf("Assay '%s' absent.", assay))
  mat <- if (inherits(obj[[assay]], "Assay5")) {
    as.matrix(LayerData(obj[[assay]], layer = layer))
  } else {
    as.matrix(GetAssayData(obj, assay = assay, slot = layer))
  }
  v <- matrixStats::rowVars(mat)
  out <- data.frame(motif = rownames(mat), variability = v, stringsAsFactors = FALSE)
  out[order(out$variability, decreasing = TRUE), , drop = FALSE]
}

chromvar_variability_plot <- function(var_df, top_n = 30) {
  stopifnot(all(c("motif","variability") %in% colnames(var_df)))
  dd <- head(var_df, top_n)
  dd$motif <- factor(dd$motif, levels = rev(dd$motif))
  ggplot(dd, aes(x = motif, y = variability)) +
    geom_col() +
    coord_flip() +
    labs(title = sprintf("chromVAR: top %d motifs par variance", nrow(dd)),
         x = "Motif", y = "Variance des déviations")
}

choose_display_reduction <- function(obj) {
  reds <- names(obj@reductions)
  prefs <- c("multimodal_umap", "RNA_umap", "atac_umap", "umap", "umap.rna", "umap.atac")
  hit <- prefs[prefs %in% reds][1]
  if (!is.na(hit)) return(hit)
  umaps <- grep("umap", reds, value = TRUE, ignore.case = TRUE)
  if (length(umaps) > 0) return(umaps[1])
  stop(sprintf("Aucune réduction UMAP trouvée. Disponibles: %s", paste(reds, collapse = ", ")))
}

do_differential_accessibility <- function(obj, peaks_assay) {
  DefaultAssay(obj) <- peaks_assay
  latent_var <- if (is.null(CONFIG$da_latent_var)) paste0("nCount_", peaks_assay) else CONFIG$da_latent_var
  if (!latent_var %in% colnames(obj@meta.data)) latent_var <- NULL
  
  idents <- levels(Idents(obj))
  if (length(idents) < 2) stop("DA exige ≥2 clusters dans Idents(obj).")
  
  res_list <- lapply(idents, function(id) {
    tbl <- FindMarkers(
      object            = obj,
      ident.1           = id,
      test.use          = CONFIG$da_test,
      latent.vars       = latent_var,
      min.cells.group   = 1,
      min.cells.feature = 1,
      min.pct           = CONFIG$da_min_pct,   # 0.05
      logfc.threshold   = 0,
      only.pos          = FALSE
    )
    tbl$cluster <- id
    tbl$feature <- rownames(tbl)
    tbl
  })
  da <- do.call(rbind, res_list)
  rownames(da) <- NULL
  da
}

do_differential_expression <- function(obj, rna_assay) {
  DefaultAssay(obj) <- rna_assay
  idents <- levels(Idents(obj))
  if (length(idents) < 2) stop("DE exige ≥2 clusters dans Idents(obj).")
  
  res_list <- lapply(idents, function(id) {
    tbl <- FindMarkers(
      object            = obj,
      ident.1           = id,
      min.cells.group   = 1,
      min.cells.feature = 1,
      min.pct           = CONFIG$de_min_pct,  # 0
      logfc.threshold   = 0,
      only.pos          = FALSE
    )
    tbl$cluster <- id
    tbl$gene <- rownames(tbl)
    tbl
  })
  de <- do.call(rbind, res_list)
  rownames(de) <- NULL
  de
}

# ------------------------------ LOAD -----------------------------------------
if (!file.exists(CONFIG$integrated_rds)) stop("Fichier RDS introuvable.")
obj <- readRDS(CONFIG$integrated_rds)
obj <- ensure_groupings(obj)

message("Assays: ", paste(Assays(obj), collapse = ", "))
peaks_assay <- choose_peaks_assay(obj)
message("Assay chromatin: '", peaks_assay, "'")
if (!(CONFIG$rna_assay %in% Assays(obj))) stop(sprintf("Assay RNA '%s' absent.", CONFIG$rna_assay))

genome <- load_bsgenome(CONFIG$bsgenome_package, seqstyle = CONFIG$seqstyle)

# ------------------------------ BLOCK 1: LinkPeaks ---------------------------
BLOCK1_OK <- FALSE
try({
  obj <- run_linkpeaks_once(obj, genome = genome, peaks_assay = peaks_assay, cfg = CONFIG$linkpeaks)
  
  link_tbl <- tryCatch(as.data.frame(Links(obj[[peaks_assay]])), error = function(e) NULL)
  if (!is.null(link_tbl) && nrow(link_tbl) > 0) {
    p_links_score <- ggplot(link_tbl, aes(x = score)) +
      geom_histogram(bins = 50) +
      labs(title = "Peak→gene link scores", x = "score", y = "count")
    p_links_pval <- if ("pvalue" %in% colnames(link_tbl)) {
      ggplot(link_tbl, aes(x = -log10(pvalue))) +
        geom_histogram(bins = 50) +
        labs(title = "Peak→gene link significance", x = "-log10(pvalue)", y = "count")
    } else NULL
    # print(p_links_score); if (!is.null(p_links_pval)) print(p_links_pval)
    
    if (isTRUE(CONFIG$do_export)) {
      dir.create(CONFIG$out_dir, showWarnings = FALSE, recursive = TRUE)
      utils::write.csv(link_tbl, file = file.path(CONFIG$out_dir, "links_table.csv"), row.names = FALSE)
      # ggsave(file.path(CONFIG$out_dir, "links_score_hist.png"), p_links_score, width = 6, height = 4, dpi = 300)
      # if (!is.null(p_links_pval)) ggsave(file.path(CONFIG$out_dir, "links_pval_hist.png"), p_links_pval, 6, 4, dpi = 300)
    }
    message(sprintf("LinkPeaks: %d liens disponibles.", nrow(link_tbl)))
  } else {
    message("LinkPeaks: aucun lien à visualiser.")
  }
  
  BLOCK1_OK <- TRUE
}, silent = FALSE)
message(if (BLOCK1_OK) "✓ BLOCK 1 (LinkPeaks) completed." else "✗ BLOCK 1 (LinkPeaks) error.")

# ============================================================
# Export LinkPeaks table to filtered_linkpeaks_all.csv
# ============================================================
message("[Export] Writing filtered_linkpeaks_all.csv ...")

# Select the correct assay containing the peaks
if ("peaks" %in% names(obj@assays)) {
  assay_name <- "peaks"
} else if ("ATAC" %in% names(obj@assays)) {
  assay_name <- "ATAC"
} else {
  stop("No 'peaks' or 'ATAC' assay found in object.")
}

# Extract Links
lk <- tryCatch({
  as.data.frame(Links(obj[[assay_name]]))
}, error = function(e) {
  stop("Cannot extract Links from assay '", assay_name, "'. Did you run LinkPeaks() earlier?")
})

# Check required columns
if (!all(c("gene", "peak") %in% colnames(lk))) {
  stop("LinkPeaks object does not contain expected columns 'gene' and 'peak'.")
}

# Split peak coordinates into chromosome/start/end
parts <- do.call(rbind, strsplit(lk$peak, "-"))
colnames(parts) <- c("Chrx", "start", "end")

# Combine and save
out <- data.frame(
  gene  = lk$gene,
  Chrx  = parts[, "Chrx"],
  start = as.integer(parts[, "start"]),
  end   = as.integer(parts[, "end"]),
  peak  = lk$peak,
  stringsAsFactors = FALSE
)

write.csv(out, "filtered_linkpeaks_all.csv", row.names = FALSE)
message("[Export] File written: ", normalizePath("filtered_linkpeaks_all.csv"))


# ------------------------------ BLOCK 2: Motifs / chromVAR -------------------
BLOCK2_OK <- FALSE
try({
  obj <- run_chromvar_once(obj, genome = genome, peaks_assay = peaks_assay)
  
  p_chromvar_var  <- NULL
  p_chromvar_umap <- NULL
  
  if ("chromvar" %in% Assays(obj)) {
    var_df <- chromvar_variability_df(obj, assay = "chromvar", layer = "data")
    p_chromvar_var <- chromvar_variability_plot(var_df, top_n = 30)
    print(p_chromvar_var)
    
    top_motif <- intersect(var_df$motif, rownames(obj[["chromvar"]]))[1]
    if (!is.na(top_motif)) {
      red <- choose_display_reduction(obj)
      old_assay <- DefaultAssay(obj); DefaultAssay(obj) <- "chromvar"
      p_chromvar_umap <- FeaturePlot(
        obj,
        reduction  = red,
        features   = top_motif,
        order      = TRUE,
        min.cutoff = "q10",
        max.cutoff = "q90"
      ) + ggplot2::ggtitle(paste("chromVAR deviation:", top_motif))
      print(p_chromvar_umap)
      DefaultAssay(obj) <- old_assay
    } else {
      message("Aucun motif top commun nom_table/assay chromVAR → skip UMAP.")
    }
  } else {
    message("Assay 'chromvar' absent → pas de variabilité/UMAP.")
  }
  
  BLOCK2_OK <- TRUE
}, silent = FALSE)
message(if (BLOCK2_OK) "✓ BLOCK 2 (Motifs/chromVAR) completed." else "✗ BLOCK 2 (Motifs/chromVAR) error.")

# -------- BLOCK 2c: chromVAR enrichments (cluster vs rest), min.pct=0.05 ----
BLOCK2c_OK <- FALSE
try({
  stopifnot("chromvar" %in% Assays(obj))
  clust_levels <- levels(Idents(obj))
  target_clusters <- clust_levels[grepl("^Cluster_", clust_levels)]
  
  latent_var <- if ("nCount_peaks" %in% colnames(obj@meta.data)) "nCount_peaks" else NULL
  DefaultAssay(obj) <- "chromvar"
  
  if (isTRUE(CONFIG$do_export)) dir.create(CONFIG$out_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (cl in target_clusters) {
    tb <- Seurat::FindMarkers(
      object            = obj,
      ident.1           = cl,
      only.pos          = FALSE,
      test.use          = "LR",
      min.pct           = CONFIG$chromvar_min_pct,  # 0.05
      latent.vars       = latent_var,
      min.cells.group   = 1,
      min.cells.feature = 1,
      logfc.threshold   = 0
    )
    tb$motif <- rownames(tb); rownames(tb) <- NULL
    
    if (isTRUE(CONFIG$do_export)) {
      utils::write.csv(tb, file = file.path(CONFIG$out_dir, sprintf("chromvar_subset_%s.csv", sub("^Cluster_", "", cl))), row.names = FALSE)
    }
  }
  
  BLOCK2c_OK <- TRUE
}, silent = FALSE)
if (BLOCK2c_OK) message("✓ BLOCK 2c (chromVAR enrichments) completed.")

# ------------------------------ BLOCK 3: DA by cluster -----------------------
BLOCK3_OK <- FALSE
da_cluster <- NULL
try({
  if (CONFIG$do_da_by_cluster) {
    DefaultAssay(obj) <- peaks_assay
    da_cluster <- do_differential_accessibility(obj, peaks_assay = peaks_assay)
    
    message(sprintf("DA: %d lignes, %d clusters.", nrow(da_cluster), length(unique(da_cluster$cluster))))
    
    example_cluster <- unique(da_cluster$cluster)[1]
    da_one <- subset(da_cluster, cluster == example_cluster)
    if (!"p_val_adj" %in% colnames(da_one) && "p_val" %in% colnames(da_one)) {
      da_one$p_val_adj <- p.adjust(da_one$p_val, method = "BH")
    }
    p_da <- ggplot(da_one, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
      geom_point(alpha = 0.5) +
      labs(title = paste("DA (cluster vs rest):", example_cluster),
           x = "avg_log2FC", y = "-log10(adj p)")
    # print(p_da)
    
    if (isTRUE(CONFIG$do_export)) {
      dir.create(CONFIG$out_dir, showWarnings = FALSE, recursive = TRUE)
      utils::write.csv(da_cluster, file = file.path(CONFIG$out_dir, "da_by_cluster.csv"), row.names = FALSE)
      # ggsave(file.path(CONFIG$out_dir, paste0("DA_volcano_", example_cluster, ".png")), p_da, 6, 5, dpi = 300)
    }
  } else {
    message("DA disabled.")
  }
  BLOCK3_OK <- TRUE
}, silent = FALSE)
message(if (BLOCK3_OK) "✓ BLOCK 3 (DA by cluster) completed." else "✗ BLOCK 3 (DA by cluster) error.")

# ------------------------------ BLOCK 4: DE by cluster -----------------------
BLOCK4_OK <- FALSE
de_cluster <- NULL
try({
  if (CONFIG$do_de_by_cluster) {
    DefaultAssay(obj) <- CONFIG$rna_assay
    de_cluster <- do_differential_expression(obj, rna_assay = CONFIG$rna_assay)
    
    message(sprintf("DE: %d lignes, %d clusters.", nrow(de_cluster), length(unique(de_cluster$cluster))))
    
    example_cluster_de <- unique(de_cluster$cluster)[1]
    de_one <- subset(de_cluster, cluster == example_cluster_de)
    if (!"p_val_adj" %in% colnames(de_one) && "p_val" %in% colnames(de_one)) {
      de_one$p_val_adj <- p.adjust(de_one$p_val, method = "BH")
    }
    p_de <- ggplot(de_one, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
      geom_point(alpha = 0.5) +
      labs(title = paste("DE (cluster vs rest):", example_cluster_de),
           x = "avg_log2FC", y = "-log10(adj p)")
    # print(p_de)
    
    if (isTRUE(CONFIG$do_export)) {
      dir.create(CONFIG$out_dir, showWarnings = FALSE, recursive = TRUE)
      utils::write.csv(de_cluster, file = file.path(CONFIG$out_dir, "de_by_cluster.csv"), row.names = FALSE)
      # ggsave(file.path(CONFIG$out_dir, paste0("DE_volcano_", example_cluster_de, ".png")), p_de, 6, 5, dpi = 300)
    }
  } else {
    message("DE disabled.")
  }
  BLOCK4_OK <- TRUE
}, silent = FALSE)
message(if (BLOCK4_OK) "✓ BLOCK 4 (DE by cluster) completed." else "✗ BLOCK 4 (DE by cluster) error.")

# ------------------------------ BLOCK 5: Subset (optionnel) ------------------
BLOCK5_OK <- FALSE
obj_subset <- NULL
try({
  target_label <- "Cluster_G"
  
  if (target_label %in% levels(Idents(obj))) {
    obj_subset <- subset(obj, idents = target_label)
  } else {
    md <- obj@meta.data
    found_col <- NULL
    for (cl in colnames(md)) {
      v <- md[[cl]]
      if ((is.factor(v) || is.character(v)) && target_label %in% as.character(v)) { found_col <- cl; break }
    }
    if (!is.null(found_col)) {
      keep_cells <- rownames(md)[as.character(md[[found_col]]) == target_label]
      obj_subset <- subset(obj, cells = keep_cells)
    } else {
      message(sprintf("Label '%s' introuvable → skip subset.", target_label))
    }
  }
  
  if (!is.null(obj_subset)) {
    use_assay <- choose_peaks_assay(obj_subset)
    DefaultAssay(obj_subset) <- use_assay
    obj_subset <- DietSeurat(obj_subset, assays = use_assay, dimreducs = NULL, graphs = NULL)
    message(sprintf("Subset '%s': %d cellules, assay '%s'.", target_label, ncol(obj_subset), use_assay))
    # if (isTRUE(CONFIG$do_export)) saveRDS(obj_subset, file = file.path(CONFIG$out_dir, "obj_subset_Cluster_G_peaks_only.rds"))
  }
  
  BLOCK5_OK <- TRUE
}, silent = FALSE)
message(if (BLOCK5_OK) "✓ BLOCK 5 (Subset) completed." else "✗ BLOCK 5 (Subset) error.")

# ------------------------------ SUMMARY --------------------------------------
message("\nSummary:")
message(" - BLOCK 1 (LinkPeaks):       ", ifelse(BLOCK1_OK, "OK", "ERROR"))
message(" - BLOCK 2 (Motifs/chromVAR): ", ifelse(BLOCK2_OK, "OK", "ERROR"))
message(" - BLOCK 2c (chromVAR enrich):", ifelse(BLOCK2c_OK, "OK", "ERROR"))
message(" - BLOCK 3 (DA by cluster):   ", ifelse(BLOCK3_OK, "OK", "ERROR"))
message(" - BLOCK 4 (DE by cluster):   ", ifelse(BLOCK4_OK, "OK", "ERROR"))
message(" - BLOCK 5 (Subset):          ", ifelse(BLOCK5_OK, "OK", "ERROR"))

# Sauvegarde finale optionnelle
# if (isTRUE(CONFIG$do_export)) saveRDS(obj, file = file.path(CONFIG$out_dir, "integrated_links_motifs_da_de.rds"))

