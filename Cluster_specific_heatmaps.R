# ------------------------------------------------------------
# 03_visualization_motifs_heatmaps.R
# Seurat/Signac multiome viz + motif summaries + heatmaps (RNA, GA, DAR, TFBS)
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(grid)
  library(viridisLite)
  library(TFBSTools)
  library(JASPAR2022)
})

set.seed(123)

# -------------------------------
# CONFIG
# -------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

CONFIG <- list(
  # The script will try these candidates in order and load the first existing file.
  seurat_rds_candidates = c(
    Sys.getenv("SEURAT_RDS", unset = NA),
    "./integrated_links_motifs_da_de.rds",       # from the new pipeline (optional)
    "./integrated_multiome_links_motifs.rds",    # legacy name (optional)
    "./integrated.rds",                          # older object (optional)
    "C:/Users/Charles/Desktop/Charles/SnMultiome/Data/PGintegrated73.RN7.links.rds" # your path (optional)
  ),
  
  out_dir            = Sys.getenv("OUT_DIR", unset = "./viz_outputs"),
  atac_umap_name     = "ATAC_umap_by_sample",
  
  # Small feature UMAPs (RNA)
  feature_genes      = c("Mettl15", "Grem1"),
  
  # Cluster order to display in heatmaps
  desired_clusters   = c("EC","Pe","Le","FSC","M","C","S","L","T","G"),
  
  # Sampling for per-cell heatmap
  per_cluster_cells  = 100,
  
  # Where per-cluster chromVAR CSVs (chromvar_subset_*.csv) live
  # If you exported them in the previous script, point here (often the same as that script's out_dir).
  chromvar_dir       = Sys.getenv("CHROMVAR_DIR", unset = "./data_processed"),
  chromvar_pattern   = "^chromvar_subset_.*\\.csv$",
  
  # Aggregated matrices for DAR/TFBS heatmaps (optional). If not present, those sections are skipped.
  fusion_peaks_csv   = Sys.getenv("FUSION_PEAKS", unset = "./fusion_peaks.csv"),
  fusion_tfbs_csv    = Sys.getenv("FUSION_TFBS",  unset = "./ChromVar/fusion_chromvar.csv"),
  
  # Thresholds for filtering aggregated matrices
  dar_padj_threshold  = 1e-7,
  tfbs_padj_threshold = 1e-3,
  
  # Rows to keep by variance for big matrices
  heatmap_top_var_n  = 5000
)

# -------------------------------
# IO helpers
# -------------------------------
.dir_create <- function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE)
.save_plot <- function(p, filename, width=14, height=10, dpi=600) {
  ggsave(filename = file.path(CONFIG$out_dir, paste0(filename, ".svg")),  plot = p, width=width, height=height)
  ggsave(filename = file.path(CONFIG$out_dir, paste0(filename, ".tiff")), plot = p, width=width, height=height,
         dpi=dpi, device="tiff", compression="lzw", units="in")
}

# Seurat v4/v5 matrix getter
get_matrix <- function(obj, assay, layer_or_slot = c("data","counts")) {
  layer_or_slot <- match.arg(layer_or_slot)
  if ("LayerData" %in% getNamespaceExports("Seurat")) {
    return(LayerData(object = obj[[assay]], layer = layer_or_slot))
  } else {
    return(GetAssayData(object = obj, assay = assay, slot = layer_or_slot))
  }
}

# UMAP chooser for display (robust to naming)
choose_display_reduction <- function(obj) {
  reds <- names(obj@reductions)
  prefs <- c("multimodal_umap","RNA_umap","atac_umap","umap","umap.rna","umap.atac")
  hit <- prefs[prefs %in% reds][1]
  if (!is.na(hit)) return(hit)
  umaps <- grep("umap", reds, value = TRUE, ignore.case = TRUE)
  if (length(umaps) > 0) return(umaps[1])
  stop(sprintf("No UMAP reduction found. Available: %s", paste(reds, collapse = ", ")))
}

# Simple log2 enrichment helper
log2_enrichment <- function(x_in, x_out, pseudo=1) log2((x_in + pseudo)/(x_out + pseudo))

# Fixed color mapping for clusters
cluster_colors <- c(
  S   = "#F8766D",
  L   = "#DB8E00",
  FSC = "#AEA200",
  C   = "#00BD5C",
  Le  = "#00C1A7",
  G   = "#00BADE",
  T   = "#00A6FF",
  M   = "#B385FF",
  EC  = "#EF67EB",
  Pe  = "#FF63B6"
)

# -------------------------------
# 1) Load Seurat object (try multiple candidates)
# -------------------------------
.dir_create(CONFIG$out_dir)

load_first_existing <- function(paths) {
  for (p in paths) {
    if (is.na(p)) next
    if (file.exists(p)) {
      message("Loading Seurat object: ", p)
      return(readRDS(p))
    }
  }
  stop("No Seurat RDS found among candidates. Please set SEURAT_RDS or update CONFIG$seurat_rds_candidates.")
}

obj <- load_first_existing(CONFIG$seurat_rds_candidates)

# Minimal metadata normalization:
# - cluster_id: prefer 'cell_Type' (same convention as the integration pipeline), else Idents
# - orig.ident/dataset for sample label
if (!"cluster_id" %in% colnames(obj@meta.data)) {
  if ("cell_Type" %in% colnames(obj@meta.data)) {
    obj$cluster_id <- as.character(obj$cell_Type)
  } else {
    obj$cluster_id <- as.character(Idents(obj))
  }
}
if (!"orig.ident" %in% colnames(obj@meta.data)) {
  obj$orig.ident <- obj$dataset %||% obj$orig.ident %||% "sample"
}
obj$cluster_clean <- gsub("^Cluster_", "", obj$cluster_id)

# -------------------------------
# 2) UMAPs: ATAC by sample + RNA Feature UMAPs
# -------------------------------
try({
  if ("atac_umap" %in% names(obj@reductions)) {
    p_atac <- DimPlot(obj, reduction = "atac_umap", group.by = "orig.ident", label = FALSE, pt.size = 0.05) +
      ggtitle("ATAC UMAP by sample") + theme_classic(base_size = 12)
    .save_plot(p_atac, CONFIG$atac_umap_name, width=14, height=10)
  }
}, silent = TRUE)

try({
  DefaultAssay(obj) <- "SoupXRNA"
  red <- tryCatch(choose_display_reduction(obj), error = function(e) names(obj@reductions)[1])
  for (g in CONFIG$feature_genes) {
    if (g %in% rownames(get_matrix(obj, "SoupXRNA", "data"))) {
      p <- FeaturePlot(obj, reduction = red, features = g, max.cutoff = "q90", pt.size = 0.05, order = TRUE) +
        ggtitle(sprintf("RNA feature: %s", g))
      .save_plot(p, paste0("UMAP_", g), width=6, height=4)
    }
  }
}, silent = TRUE)

# -------------------------------
# 3) Motif family summary per cluster from chromVAR CSVs (optional)
# -------------------------------
try({
  opts <- list(collection = "CORE", tax_group = "vertebrates")
  motifSet <- TFBSTools::getMatrixSet(JASPAR2022, opts)
  jaspar_lookup <- data.frame(
    MatrixID = vapply(motifSet, function(m) m@ID, FUN.VALUE = character(1)),
    Family   = vapply(motifSet, function(m) {
      fam <- m@tags$family
      if (is.null(fam) || length(fam) == 0) NA_character_ else fam[[1]]
    }, FUN.VALUE = character(1)),
    stringsAsFactors = FALSE
  )
  
  files <- list.files(path = CONFIG$chromvar_dir, pattern = CONFIG$chromvar_pattern, full.names = TRUE)
  if (length(files) > 0) {
    sink(file.path(CONFIG$out_dir, "chromvar_top3_families.txt"))
    for (f in files) {
      df <- read.csv(f, stringsAsFactors = FALSE)
      # Support both column names: 'X' (motif id) or 'motif'
      motif_col <- if ("X" %in% names(df)) "X" else if ("motif" %in% names(df)) "motif" else NULL
      if (is.null(motif_col)) next
      # Rank by adjusted p then effect size
      ord <- order(df$p_val_adj, -df$avg_log2FC)
      df <- df[ord, ]
      top10 <- head(df, 10)
      merged <- merge(top10, jaspar_lookup, by.x = motif_col, by.y = "MatrixID", all.x = TRUE)
      family_summary <- merged %>%
        group_by(Family) %>%
        summarise(min_p_val_adj = min(p_val_adj, na.rm = TRUE), .groups = "drop") %>%
        arrange(min_p_val_adj)
      cat("File:", basename(f), "\n")
      print(head(as.data.frame(family_summary), 3)); cat("\n---\n\n")
    }
    sink()
  }
}, silent = TRUE)

# Also capture first motif / avg_log2FC / first non-zero pval per file (optional)
try({
  files <- list.files(path = CONFIG$chromvar_dir, pattern = CONFIG$chromvar_pattern, full.names = TRUE)
  if (length(files) > 0) {
    out_rows <- lapply(files, function(f) {
      df <- read.csv(f, stringsAsFactors = FALSE)
      motif_col <- if ("X" %in% names(df)) "X" else if ("motif" %in% names(df)) "motif" else NULL
      if (is.null(motif_col)) return(NULL)
      ord <- order(df$p_val_adj, -df$avg_log2FC)
      df <- df[ord, ]
      motif <- df[[motif_col]][1]
      avg_log2FC_val <- df$avg_log2FC[1]
      nz <- which(df$p_val_adj != 0)
      first_nz <- if (length(nz) > 0) df$p_val_adj[nz[1]] else NA
      data.frame(file = basename(f), motif = motif, avg_log2FC = avg_log2FC_val, first_nonzero_padj = first_nz)
    })
    out_rows <- out_rows[!vapply(out_rows, is.null, logical(1))]
    if (length(out_rows) > 0) {
      write.csv(bind_rows(out_rows), file.path(CONFIG$out_dir, "chromvar_first_motif_per_cluster.csv"), row.names = FALSE)
    }
  }
}, silent = TRUE)

# -------------------------------
# 4) Per-cell RNA “marker panel” heatmap (sampled cells per cluster)
# -------------------------------
try({
  DefaultAssay(obj) <- "SoupXRNA"
  
  # Ensure desired clusters exist
  clusters_present <- intersect(CONFIG$desired_clusters, unique(obj$cluster_clean))
  stopifnot(length(clusters_present) > 0)
  
  # Sample cells per cluster (safe)
  selected_cells <- unlist(lapply(unique(obj$cluster_id), function(clust) {
    cells <- rownames(obj@meta.data)[obj$cluster_id == clust]
    if (length(cells) == 0) return(character(0))
    sample(cells, size = min(CONFIG$per_cluster_cells, length(cells)))
  }))
  
  # Panel of genes (customize as you wish)
  genes_of_interest <- c(
    "Emcn","Plvap","Pde5a","Pdgfrb","Tyrobp","Arhgap15",
    "Adamts9","Rfx4","S100b","Pomc","Neurod1",
    "Pax7","Pou1f1","Prl","S100g","Gh1","Tshb","Cga",
    "Lhb","Fshb","Nr5a1","Gnrhr"
  )
  
  rna_counts <- get_matrix(obj, "SoupXRNA", "counts")
  genes_present <- intersect(genes_of_interest, rownames(rna_counts))
  if (length(genes_present) > 0 && length(selected_cells) > 0) {
    full_counts <- rna_counts[genes_present, , drop = FALSE]
    pseudo <- 1
    log_fc <- sapply(selected_cells, function(cell) {
      cur_cl <- obj@meta.data[cell, "cluster_id"]
      cells_in  <- rownames(obj@meta.data)[obj$cluster_id == cur_cl]
      cells_out <- setdiff(colnames(full_counts), cells_in)
      cell_counts  <- full_counts[, cell]
      mean_non     <- rowMeans(full_counts[, cells_out, drop = FALSE])
      log2((cell_counts + pseudo) / (mean_non + pseudo))
    })
    log_fc <- as.matrix(log_fc)
    
    # Order columns by desired cluster sequence
    cln <- gsub("^Cluster_", "", obj$cluster_id[selected_cells])
    cl_factor <- factor(cln, levels = CONFIG$desired_clusters)
    ord <- order(cl_factor)
    log_fc_ord <- log_fc[, ord, drop = FALSE]
    clusters_ord <- cl_factor[ord]
    
    # Annotation and colors
    annotation_col <- data.frame(Cluster = clusters_ord)
    rownames(annotation_col) <- colnames(log_fc_ord)
    annotation_colors <- list(Cluster = cluster_colors[CONFIG$desired_clusters])
    
    heatmap_colors <- colorRampPalette(c("beige", "brown2"))(100)
    my_breaks <- seq(0, 3, length.out = 101)
    
    hm1 <- pheatmap(
      log_fc_ord, scale = "row", cluster_rows = FALSE, cluster_cols = FALSE,
      show_colnames = FALSE, annotation_col = annotation_col,
      annotation_colors = annotation_colors, color = heatmap_colors,
      breaks = my_breaks, fontsize_row = 10,
      main = "Per-cell log2 enrichment (RNA markers)"
    )
    
    svg(file.path(CONFIG$out_dir, "heatmap_markers_cells.svg"), width = 10, height = 6)
    grid.newpage(); grid.draw(hm1$gtable); dev.off()
    
    tiff(file.path(CONFIG$out_dir, "heatmap_markers_cells.tiff"), units = "in", width = 10, height = 6, res = 600)
    grid.newpage(); grid.draw(hm1$gtable); dev.off()
  }
}, silent = TRUE)

# -------------------------------
# 5) Top-10 RNA markers per cluster heatmap (cluster means, log2 enrich)
# -------------------------------
try({
  DefaultAssay(obj) <- "SoupXRNA"
  
  markers <- FindAllMarkers(object = obj, assay = "SoupXRNA", only.pos = TRUE,
                            logfc.threshold = 0.25, min.pct = 0.1)
  
  # Select top-10 per cluster (filter out LOC* if present)
  exclude_pattern <- "^LOC\\d"
  selected_mapping <- list()
  for (cl in CONFIG$desired_clusters) {
    cl_label <- paste0("Cluster_", cl)
    cluster_markers <- markers %>%
      filter(cluster == cl_label, !grepl(exclude_pattern, gene)) %>%
      arrange(desc(avg_log2FC))
    selected_mapping[[cl]] <- head(cluster_markers$gene, 10)
  }
  selected_genes <- unique(unlist(selected_mapping, use.names = FALSE))
  if (length(selected_genes) > 0) {
    # Per-cluster means
    rna_counts <- get_matrix(obj, "SoupXRNA", "counts")
    mat_counts <- rna_counts[selected_genes, , drop = FALSE]
    cells_by_cluster <- split(rownames(obj@meta.data), obj$cluster_clean)
    
    mean_expr <- matrix(0, nrow = length(selected_genes), ncol = length(CONFIG$desired_clusters),
                        dimnames = list(selected_genes, CONFIG$desired_clusters))
    for (cl in CONFIG$desired_clusters) {
      if (!cl %in% names(cells_by_cluster)) next
      cs <- cells_by_cluster[[cl]]
      mean_expr[, cl] <- rowMeans(mat_counts[, cs, drop = FALSE])
    }
    
    # Log2 enrichment vs rest
    pseudo <- 1
    heatmap_matrix <- matrix(0, nrow = length(selected_genes), ncol = length(CONFIG$desired_clusters),
                             dimnames = list(selected_genes, CONFIG$desired_clusters))
    for (g in selected_genes) {
      for (cl in CONFIG$desired_clusters) {
        in_cl  <- mean_expr[g, cl]
        out_cl <- mean(colMeans(as.matrix(mean_expr[g, setdiff(CONFIG$desired_clusters, cl), drop=TRUE])))
        heatmap_matrix[g, cl] <- log2_enrichment(in_cl, out_cl, pseudo)
      }
    }
    
    final_colors <- colorRampPalette(c("blue","white","red"))(100)
    my_breaks2 <- seq(-8, 8, length.out = 101)
    
    hm2 <- pheatmap(
      heatmap_matrix, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE,
      annotation_col = data.frame(Cluster = factor(CONFIG$desired_clusters, levels = CONFIG$desired_clusters),
                                  row.names = CONFIG$desired_clusters),
      annotation_colors = list(Cluster = cluster_colors[CONFIG$desired_clusters]),
      color = final_colors, breaks = my_breaks2,
      main = "Top10 RNA markers per cluster (log2 enrich)", fontsize_row = 9, cellwidth = 24
    )
    
    svg(file.path(CONFIG$out_dir, "heatmap_markers_top10.svg"), width = 13, height = 20)
    grid.newpage(); grid.draw(hm2$gtable); dev.off()
    
    tiff(file.path(CONFIG$out_dir, "heatmap_markers_top10.tiff"), units = "in", width = 13, height = 20, res = 600)
    grid.newpage(); grid.draw(hm2$gtable); dev.off()
  }
}, silent = TRUE)

# -------------------------------
# 6) GeneActivity enrichment heatmap (markers, cluster vs rest)
# -------------------------------

try({
  # --- Parameters you may edit ---
  desired_cluster_order <- c("EC","Pe","Le","FSC","M","C","S","L","T","G")
  genes_of_interest <- c(
    "Emcn","Plvap","Pde5a","Pdgfrb","Tyrobp",
    "Adamts9","Rfx4","S100b","Prl","S100g",
    "Pou1f1","Gh1","Tshb","Cga","Pomc","Neurod1",
    "Pax7","Lhb","Fshb","Nr5a1","Gnrhr"
  )
  out_svg  <- file.path(CONFIG$out_dir, "heatmap_geneactivity.svg")
  out_tiff <- file.path(CONFIG$out_dir, "heatmap_geneactivity.tiff")
  
  # --- Sanity: assay present ---
  if (!"GeneActivity" %in% names(obj@assays)) stop("GeneActivity assay not found.")
  
  # --- Robust matrix getter (Seurat v4/v5) ---
  ga_counts <- get_matrix(obj, assay = "GeneActivity", layer_or_slot = "counts")
  all_genes <- rownames(ga_counts)
  
  # --- Keep only available markers; warn on missing ---
  genes_present <- intersect(genes_of_interest, all_genes)
  if (length(genes_present) == 0L) stop("None of the marker genes are present in GeneActivity.")
  missing_genes <- setdiff(genes_of_interest, all_genes)
  if (length(missing_genes) > 0L) {
    message("[GA] Missing in GeneActivity: ", paste(missing_genes, collapse = ", "))
  }
  
  # --- Cluster labels: ensure 'cluster_clean' exists (no 'Cluster_' prefix) ---
  if (!"cluster_id" %in% colnames(obj@meta.data)) {
    obj$cluster_id <- if ("cell_Type" %in% colnames(obj@meta.data)) obj$cell_Type else Idents(obj)
  }
  obj$cluster_clean <- gsub("^Cluster_", "", as.character(obj$cluster_id))
  
  # --- Build enrichment matrix (log2 in/out with pseudo=1) ---
  pseudo <- 1
  enrichment_matrix <- matrix(
    NA_real_,
    nrow = length(genes_present),
    ncol = length(desired_cluster_order),
    dimnames = list(genes_present, desired_cluster_order)
  )
  
  cells_by_cl <- split(rownames(obj@meta.data), obj$cluster_clean)
  for (cl in desired_cluster_order) {
    if (!cl %in% names(cells_by_cl)) next
    cells_in  <- cells_by_cl[[cl]]
    cells_out <- setdiff(colnames(ga_counts), cells_in)
    
    mean_in  <- rowMeans(ga_counts[genes_present, cells_in,  drop = FALSE])
    mean_out <- rowMeans(ga_counts[genes_present, cells_out, drop = FALSE])
    
    enrichment_matrix[, cl] <- log2((mean_in + pseudo) / (mean_out + pseudo))
  }
  
  # --- Keep genes in desired display order (only those present) ---
  desired_gene_order <- c(
    "Emcn","Plvap","Pde5a","Pdgfrb","Tyrobp",
    "Adamts9","Rfx4","S100b","Pomc","Neurod1",
    "Pax7","Pou1f1","Prl","S100g","Gh1","Tshb","Cga",
    "Lhb","Fshb","Nr5a1","Gnrhr"
  )
  desired_gene_order <- intersect(desired_gene_order, rownames(enrichment_matrix))
  enrichment_matrix  <- enrichment_matrix[desired_gene_order, , drop = FALSE]
  
  # --- Column annotation/colors (fixed mapping) ---
  cluster_color_mapping <- c(
    S   = "#F8766D",
    L   = "#DB8E00",
    FSC = "#AEA200",
    C   = "#00BD5C",
    Le  = "#00C1A7",
    G   = "#00BADE",
    T   = "#00A6FF",
    M   = "#B385FF",
    EC  = "#EF67EB",
    Pe  = "#FF63B6"
  )
  annotation_col    <- data.frame(Cluster = colnames(enrichment_matrix),
                                  row.names = colnames(enrichment_matrix))
  annotation_colors <- list(Cluster = cluster_color_mapping[colnames(enrichment_matrix)])
  
  # --- Display palette/breaks (same as your spec) ---
  heatmap_colors_ga <- colorRampPalette(c("white", "#1aa14e"))(100)
  my_breaks_ga      <- seq(0, 8, length.out = 101)
  
  # --- Draw (no clustering; fixed column order) ---
  hm_ga <- pheatmap::pheatmap(
    enrichment_matrix,
    scale = "row",                 # row-wise z-scoring for display
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    color = heatmap_colors_ga,
    breaks = my_breaks_ga,
    fontsize_row = 12,
    main = "GeneActivity enrichment (log2 vs rest)"
  )
  print(hm_ga)
  
  # --- Save SVG and TIFF ---
  svg(out_svg, width = 10, height = 6)
  grid::grid.newpage(); grid::grid.draw(hm_ga$gtable); dev.off()
  
  tiff(out_tiff, units = "in", width = 10, height = 6, res = 600, compression = "lzw")
  grid::grid.newpage(); grid::grid.draw(hm_ga$gtable); dev.off()
}, silent = FALSE)

# -------------------------------
# 7) DAR heatmap from fusion_peaks.csv 
# -------------------------------
try({
  message("[ DAR ] Reading: ", CONFIG$fusion_peaks_csv)
  stopifnot(file.exists(CONFIG$fusion_peaks_csv))
  data <- read.csv(CONFIG$fusion_peaks_csv, check.names = FALSE)
  
  # --- PARAMETERS (identical to original) ---
  pval_threshold <- 1e-7
  top_n_rows     <- 500
  
  # --- Pair avg_log2FC_* with p_val_adj_* by suffix (robust to column order) ---
  cols_avg  <- grep("^avg_log2FC_", colnames(data), value = TRUE)
  cols_pval <- grep("^p_val_adj_",  colnames(data), value = TRUE)
  suf_avg   <- sub("^avg_log2FC_", "", cols_avg)
  suf_pval  <- sub("^p_val_adj_",  "", cols_pval)
  common    <- intersect(suf_avg, suf_pval)
  if (length(common) == 0L) stop("[ DAR ] No matched avg_log2FC_*/p_val_adj_* columns.")
  
  cols_avg  <- paste0("avg_log2FC_", common)
  cols_pval <- paste0("p_val_adj_",  common)
  
  # --- Build numeric matrix and mask by adjusted p-value ---
  M <- as.matrix(data[, cols_avg,  drop = FALSE]); storage.mode(M) <- "double"
  P <- as.matrix(data[, cols_pval, drop = FALSE]); storage.mode(P) <- "double"
  M[P >= pval_threshold] <- NA_real_
  
  if (!"peak" %in% names(data)) stop("[ DAR ] Missing 'peak' column.")
  rownames(M) <- data$peak
  
  # --- Drop constant rows/cols (variance = 0 ignoring NA) ---
  row_sd <- apply(M, 1, sd, na.rm = TRUE)
  col_sd <- apply(M, 2, sd, na.rm = TRUE)
  M <- M[row_sd > 0, , drop = FALSE]
  M <- M[, col_sd > 0, drop = FALSE]
  
  # --- Column-wise z-score (mean/sd computed on non-NA) ---
  col_means <- apply(M, 2, function(x) mean(x, na.rm = TRUE))
  col_sds   <- apply(M, 2, function(x) sd(x,   na.rm = TRUE))
  col_sds[col_sds == 0 | is.na(col_sds)] <- 1
  M <- sweep(M, 2, col_means, "-")
  M <- sweep(M, 2, col_sds,   "/")
  M[is.na(M) | is.infinite(M)] <- 0
  
  # --- Keep top rows by variance (exactly 500 like the original) ---
  row_var <- apply(M, 1, var)
  keep_ix <- order(row_var, decreasing = TRUE)[seq_len(min(top_n_rows, nrow(M)))]
  H <- M[keep_ix, , drop = FALSE]
  
  # --- Show and save heatmap (no forced column order) ---
  breaks <- seq(-5, 5, length.out = 51)
  hm4 <- pheatmap::pheatmap(
    H,
    cluster_rows = TRUE,
    cluster_cols = TRUE,      # no fixed order
    show_rownames = FALSE,
    show_colnames = TRUE,
    color = colorRampPalette(c("black","white","#1aa14e"))(50),
    breaks = breaks,
    main = "DAR heatmap (column z-score; adj p < 1e-7)",
    fontsize_col = 14,
    border_color = NA
  )
  print(hm4)
  
  tiff(file.path(CONFIG$out_dir, "Heatmap_DAR_log2FC.tiff"),
       units = "in", width = 20, height = 20, res = 300)
  grid::grid.newpage(); grid::grid.draw(hm4$gtable); dev.off()
}, silent = FALSE)

# -------------------------------
# 8) TFBS (chromVAR) heatmap from fusion_chromvar.csv (same params as original)
# -------------------------------
try({
  message("[ TFBS ] Reading: ", CONFIG$fusion_tfbs_csv)
  stopifnot(file.exists(CONFIG$fusion_tfbs_csv))
  raw <- read.csv(CONFIG$fusion_tfbs_csv, check.names = FALSE, stringsAsFactors = FALSE)
  
  # Some exports have a duplicated header row; drop it if detected
  if (any(names(raw) == raw[1,1])) raw <- raw[-1, ]
  
  # --- PARAMETERS (identical to original) ---
  pval_threshold <- 1e-3
  
  # --- Coerce numerics (keep ID column 'A' if present) ---
  id_col <- if ("A" %in% names(raw)) "A" else colnames(raw)[1]
  num_cols <- setdiff(names(raw), id_col)
  suppressWarnings({
    raw[, num_cols] <- lapply(raw[, num_cols, drop = FALSE], as.numeric)
  })
  
  # --- Pair avg_log2FC_* with p_val_adj_* by suffix ---
  cols_avg  <- grep("^avg_log2FC_", names(raw), value = TRUE)
  cols_pval <- grep("^p_val_adj_",  names(raw), value = TRUE)
  suf_avg   <- sub("^avg_log2FC_", "", cols_avg)
  suf_pval  <- sub("^p_val_adj_",  "", cols_pval)
  common    <- intersect(suf_avg, suf_pval)
  if (length(common) == 0L) stop("[ TFBS ] No matched avg_log2FC_*/p_val_adj_* columns.")
  
  cols_avg  <- paste0("avg_log2FC_", common)
  cols_pval <- paste0("p_val_adj_",  common)
  
  # --- Build matrix, mask by padj, column-wise z-score ---
  M <- as.matrix(raw[, cols_avg,  drop = FALSE]); storage.mode(M) <- "double"
  P <- as.matrix(raw[, cols_pval, drop = FALSE]); storage.mode(P) <- "double"
  M[P >= pval_threshold] <- NA_real_
  rownames(M) <- raw[[id_col]]
  
  # Drop constants
  row_sd <- apply(M, 1, sd, na.rm = TRUE)
  col_sd <- apply(M, 2, sd, na.rm = TRUE)
  M <- M[row_sd > 0, , drop = FALSE]
  M <- M[, col_sd > 0, drop = FALSE]
  
  # Column z-score
  col_means <- apply(M, 2, function(x) mean(x, na.rm = TRUE))
  col_sds   <- apply(M, 2, function(x) sd(x,   na.rm = TRUE))
  col_sds[col_sds == 0 | is.na(col_sds)] <- 1
  M <- sweep(M, 2, col_means, "-")
  M <- sweep(M, 2, col_sds,   "/")
  M[is.na(M) | is.infinite(M)] <- 0
  
  # --- Show and save heatmap (no forced column order) ---
  colnames(M) <- sub("^avg_log2FC_", "", colnames(M))
  breaks <- seq(-3, 3, length.out = 51)
  hm5 <- pheatmap::pheatmap(
    M,
    cluster_rows = TRUE,
    cluster_cols = TRUE,      # no fixed order
    show_rownames = FALSE,
    show_colnames = TRUE,
    clustering_distance_rows = "manhattan",
    clustering_method = "ward.D2",
    color = colorRampPalette(c("black","white","#1aa14e"))(50),
    breaks = breaks,
    main = "TFBS (chromVAR) avg_log2FC (column z-score; adj p < 1e-3)",
    fontsize_col = 14,
    border_color = NA
  )
  print(hm5)
  
  total_w <- sum(grid::convertWidth(hm5$gtable$widths,  "in", valueOnly = TRUE))
  total_h <- sum(grid::convertHeight(hm5$gtable$heights, "in", valueOnly = TRUE))
  tiff(file.path(CONFIG$out_dir, "Heatmap_TFBS.tiff"),
       units = "in", width = max(12, total_w), height = max(12, total_h), res = 300)
  grid::grid.newpage(); grid::grid.draw(hm5$gtable); dev.off()
}, silent = FALSE)

message("✅ Completed: outputs in ", normalizePath(CONFIG$out_dir))
