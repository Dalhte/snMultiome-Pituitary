# ==============================================================================
# Subset of DEGs per cluster (clean, portable, runs after previous scripts)
# - Consumes the Seurat object `obj` already built in earlier steps
# - Works with cluster labels like "Cluster_FSC", "Cluster_G", etc.
# - Writes tidy CSVs with an explicit Gene column
# - Builds combined table, applies thresholds, and plots two heatmaps
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(pheatmap)
  library(Matrix)
})

# -----------------------------
# USER CONFIG
# -----------------------------
cluster_of_interest <- "Cluster_L"                             # e.g., "Cluster_G", "Cluster_FSC"
output_dir          <- "./deg_outputs"                           # output directory
sample_targets      <- c(D2 = 514, PM = 42, PS = 513, E = 78)    # downsampling per origin

# -----------------------------
# PRECHECKS
# -----------------------------
stopifnot(exists("obj"), inherits(obj, "Seurat"))
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!("SoupXRNA" %in% Assays(obj))) stop("Assay 'SoupXRNA' not found in `obj`.")

# Normalize identities/metadata as in previous scripts
if (!"cluster_id" %in% colnames(obj@meta.data)) {
  obj$cluster_id <- if ("cell_Type" %in% colnames(obj@meta.data)) as.character(obj$cell_Type) else as.character(Idents(obj))
}
# Ensure clean version without "Cluster_" prefix also exists (optional downstream)
obj$cluster_clean <- gsub("^Cluster_", "", obj$cluster_id)

# `orig.ident` was used for sample labels in the prior pipeline
if (!"orig.ident" %in% colnames(obj@meta.data)) {
  obj$orig.ident <- obj$dataset %||% obj$orig.ident %||% "sample"
}

# -----------------------------
# SUBSET TO CLUSTER OF INTEREST
# -----------------------------
DefaultAssay(obj) <- "SoupXRNA"

# If Idents are not the requested cluster labels, set them
if (!(cluster_of_interest %in% levels(Idents(obj)))) {
  Idents(obj) <- factor(obj$cluster_id)
}
if (!(cluster_of_interest %in% levels(Idents(obj)))) {
  stop(sprintf("Cluster '%s' not found in Idents(obj) or obj$cluster_id.", cluster_of_interest))
}

subset_Cluster_obj <- subset(obj, idents = cluster_of_interest)
if (ncol(subset_Cluster_obj) == 0L) stop("Subset is empty for cluster: ", cluster_of_interest)

# -----------------------------
# ROBUST PER-ORIGIN DOWNSAMPLING
# -----------------------------
sample_cells <- function(obj_, origin_label, n) {
  cells <- rownames(obj_@meta.data)[as.character(obj_$orig.ident) == origin_label]
  if (length(cells) == 0) return(character(0))
  sample(cells, size = min(n, length(cells)), replace = FALSE)
}

set.seed(123)
selected_cells <- unlist(mapply(
  FUN = sample_cells,
  origin_label = names(sample_targets),
  n = unname(sample_targets),
  MoreArgs = list(obj_ = subset_Cluster_obj),
  SIMPLIFY = FALSE,
  USE.NAMES = FALSE
))

# If none selected, keep all; else subset to selection
if (length(selected_cells) > 0) {
  subset_Cluster_obj <- subset(subset_Cluster_obj, cells = selected_cells)
}
if (ncol(subset_Cluster_obj) == 0L) stop("No cells remain after downsampling filter.")

# -----------------------------
# GROUP IDS FOR CONTRASTS
# -----------------------------
subset_Cluster_obj$group_id <- paste0(subset_Cluster_obj$orig.ident, "_", cluster_of_interest)
Idents(subset_Cluster_obj)  <- "group_id"

# Only define contrasts present in the data
all_groups <- levels(Idents(subset_Cluster_obj))
make_pair <- function(a, b) c(paste0(a, "_", cluster_of_interest), paste0(b, "_", cluster_of_interest))
contrast_defs <- list(
  PM_D2 = make_pair("PM", "D2"),
  PS_PM = make_pair("PS", "PM"),
  E_PS  = make_pair("E",  "PS"),
  D2_E  = make_pair("D2", "E")
)
# Keep contrasts where both groups exist
contrast_defs <- Filter(function(pair) all(pair %in% all_groups), contrast_defs)
if (length(contrast_defs) == 0L) stop("No valid group contrasts for the selected cluster.")

# -----------------------------
# FINDMARKERS PER CONTRAST
# -----------------------------
fm_args <- list(
  min.cells.group   = 1,
  min.cells.feature = 1,
  min.pct           = 0,
  logfc.threshold   = 0,
  only.pos          = FALSE
)

run_findmarkers <- function(id1, id2) {
  do.call(FindMarkers, c(list(object = subset_Cluster_obj, ident.1 = id1, ident.2 = id2), fm_args))
}

fm_list <- lapply(contrast_defs, \(ids) run_findmarkers(ids[1], ids[2]))

# -----------------------------
# WRITE CSVs WITH EXPLICIT GENE COLUMN
# -----------------------------
write_contrast_csv <- function(tab, label) {
  out <- file.path(output_dir, sprintf("Combined_%s_%s.csv", cluster_of_interest, label))
  dd  <- cbind(Gene = rownames(tab), tab)
  write.csv(dd, out, row.names = FALSE)
  invisible(out)
}
contrast_csv_paths <- mapply(write_contrast_csv, fm_list, names(fm_list), SIMPLIFY = TRUE)

# -----------------------------
# ROBUST READER FOR CONTRAST TABLES
# -----------------------------
read_and_rename <- function(path, comp_label) {
  if (!file.exists(path)) return(NULL)
  df <- tryCatch(read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE),
                 error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0L) return(NULL)
  
  # gene column candidates; fallback to first column
  gene_candidates <- c("Gene","gene","SYMBOL","symbol","feature","features","X","V1")
  gene_col <- gene_candidates[gene_candidates %in% names(df)][1]
  if (is.na(gene_col)) gene_col <- names(df)[1]
  genes <- as.character(df[[gene_col]])
  
  # reject numeric index masquerading as genes
  if (all(suppressWarnings(!is.na(as.integer(genes)))) &&
      identical(as.integer(genes), seq_len(nrow(df)))) {
    stop("First column looks like an index (no gene names) in: ", basename(path))
  }
  
  # avg_log2FC
  log2fc_col <- if ("avg_log2FC" %in% names(df)) "avg_log2FC" else {
    cand <- grep("log2FC|logFC|avg.*log2", names(df), ignore.case = TRUE, value = TRUE)
    if (length(cand) == 0) stop("No log2FC-like column in: ", basename(path))
    cand[1]
  }
  log2fc <- suppressWarnings(as.numeric(df[[log2fc_col]]))
  
  # p_val_adj (or compute from p_val)
  if ("p_val_adj" %in% names(df)) {
    padj <- suppressWarnings(as.numeric(df$p_val_adj))
  } else if ("p_val" %in% names(df)) {
    padj <- p.adjust(suppressWarnings(as.numeric(df$p_val)), method = "BH")
  } else {
    stop("No p_val_adj or p_val column in: ", basename(path))
  }
  
  out <- data.frame(
    Gene = genes,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out[[paste0("avg_log2FC_", comp_label)]] <- log2fc
  out[[paste0("p_val_adj_",  comp_label)]] <- padj
  
  # drop obvious non-genes
  out <- out[!grepl("^(LOC|NEWGENE|RGD)", out$Gene), , drop = FALSE]
  out
}

tbls_raw <- Map(read_and_rename, contrast_csv_paths, names(fm_list))
per_contrast_tables <- Filter(function(x) !is.null(x) && nrow(x) > 0, tbls_raw)
if (length(per_contrast_tables) == 0L) stop("All contrast tables are empty. Check FindMarkers outputs.")

# -----------------------------
# MERGE ALL CONTRASTS BY GENE
# -----------------------------
combined_df <- Reduce(function(x, y) dplyr::full_join(x, y, by = "Gene"), per_contrast_tables)
combined_df <- combined_df %>% dplyr::filter(!is.na(Gene)) %>% dplyr::distinct(Gene, .keep_all = TRUE)
write.csv(combined_df, file.path(output_dir, sprintf("Combined_%s.csv", cluster_of_interest)), row.names = FALSE)

# -----------------------------
# SELECT DEGS BY THRESHOLDS
# -----------------------------
avg_log2FC_columns <- grep("^avg_log2FC_", colnames(combined_df), value = TRUE)
p_val_adj_columns  <- grep("^p_val_adj_",  colnames(combined_df), value = TRUE)

fc_thr   <- 1
padj_thr <- 1e-5

condition_met <- combined_df %>%
  rowwise() %>%
  mutate(
    avg_log2FC_condition = any(abs(c_across(all_of(avg_log2FC_columns))) > fc_thr, na.rm = TRUE),
    p_val_adj_condition  = any(c_across(all_of(p_val_adj_columns)) < padj_thr,  na.rm = TRUE)
  ) %>%
  filter(avg_log2FC_condition & p_val_adj_condition) %>%
  ungroup()

write.csv(condition_met, file.path(output_dir, sprintf("DEGs_Subset_%s.csv", cluster_of_interest)), row.names = FALSE)
selected_genes <- condition_met$Gene
if (length(selected_genes) == 0L) {
  warning("No genes met the thresholds; heatmaps will be empty or minimal.")
}

# ---- Heatmap 1: log2FC of average counts per origin (within cluster) ----
avg_expr <- AverageExpression(
  subset_Cluster_obj,
  assays  = "SoupXRNA",
  slot    = "counts",
  group.by = "orig.ident"
)$SoupXRNA

# Conserver l'ordre temporel s'il existe
time_order <- intersect(c("D2","PM","PS","E"), colnames(avg_expr))
avg_expr   <- as.matrix(avg_expr[, time_order, drop = FALSE])

if (ncol(avg_expr) >= 2) {
  eps <- 1e-8
  # Matrice de sortie avec dimnames corrects dès le départ
  log2fc_matrix <- matrix(NA_real_,
                          nrow = nrow(avg_expr),
                          ncol = ncol(avg_expr),
                          dimnames = list(rownames(avg_expr), colnames(avg_expr)))
  # log2( groupe / moyenne_des_autres )
  rs <- rowSums(avg_expr)  # pour accélérer
  for (j in seq_len(ncol(avg_expr))) {
    others_mean <- (rs - avg_expr[, j, drop = FALSE]) / (ncol(avg_expr) - 1)
    log2fc_matrix[, j] <- log2((avg_expr[, j] + eps) / (others_mean + eps))
  }
  
  # Garder uniquement les gènes sélectionnés et lignes finies
  keep_rows <- intersect(selected_genes, rownames(log2fc_matrix))
  hm1_mat   <- log2fc_matrix[keep_rows, colnames(avg_expr), drop = FALSE]
  hm1_mat   <- hm1_mat[apply(hm1_mat, 1, function(x) all(is.finite(x))), , drop = FALSE]
  
  if (nrow(hm1_mat) > 1) {
    svg(file.path(output_dir, sprintf("heatmap_log2fc_counts_%s.svg", cluster_of_interest)),
        width = 10, height = 8)
    pheatmap(
      hm1_mat,
      scale         = "row",
      main          = sprintf("Log2 FC of counts vs other timepoints — %s", cluster_of_interest),
      cluster_cols  = FALSE,
      cluster_rows  = TRUE,
      show_rownames = FALSE,
      fontsize_row  = 10
    )
    dev.off()
  } else {
    message("No rows left for Heatmap 1 after filtering.")
  }
} else {
  message("Only one origin present in the subset; skipping Heatmap 1.")
}

# -----------------------------
# HEATMAP 2: avg_log2FC ACROSS COMPARISONS (PM_D2, PS_PM, E_PS)
# -----------------------------
comparisons_order <- c("avg_log2FC_PM_D2", "avg_log2FC_PS_PM", "avg_log2FC_E_PS")
existing_cols     <- intersect(comparisons_order, colnames(combined_df))

hm2_mat <- NULL
if (length(existing_cols) > 0L && length(selected_genes) > 0L) {
  idx <- match(selected_genes, combined_df$Gene)
  idx <- idx[!is.na(idx)]
  if (length(idx) > 0L) {
    hm2_mat <- as.matrix(combined_df[idx, existing_cols, drop = FALSE])
    rownames(hm2_mat) <- combined_df$Gene[idx]
    hm2_mat <- hm2_mat[apply(hm2_mat, 1, function(x) all(is.finite(as.numeric(x)))), , drop = FALSE]
  }
}

if (!is.null(hm2_mat) && nrow(hm2_mat) > 1) {
  svg(file.path(output_dir, sprintf("heatmap_comparisons_%s.svg", cluster_of_interest)), width = 8, height = 10)
  pheatmap(
    hm2_mat,
    scale         = "row",
    main          = sprintf("Avg_log2FC across comparisons — %s", cluster_of_interest),
    cluster_cols  = FALSE,
    cluster_rows  = TRUE,
    show_rownames = FALSE,
    fontsize_row  = 8
  )
  dev.off()
}

message("Done. Outputs in: ", normalizePath(output_dir))
