# =============================================================================
# Heatmaps & multi-omic integration (RNA expression, RNA velocity, ATAC)
# – ATAC robust build + spline scan + BIC segmentation integrated
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(pheatmap)
  library(viridis)
  library(segmented)
  library(Seurat)              # Seurat object / assays
  library(Signac)              # ChromatinAssay, peak ops
  library(SummarizedExperiment)
  library(data.table)
  library(ggplot2)
  library(scMEGA)
})

# --------------------------- User config --------------------------------------
paths <- list(
  integrated_rds   = "integrated.rds",
  vel_csv          = "velocity_data_all_regtime.csv",
  expr_csv         = "data_regtime.csv",
  merged_csv       = "pseudotime_regulatory/merged_regulatory_time.csv",
  links_activ_csv  = "Peak_Gene_TAD_links_activator_FDR001_empirical_spline.csv",
  links_sil_csv    = "Peak_Gene_TAD_links_silencer_FDR001_empirical_spline.csv",
  out_dir          = "outputs/heatmaps"
)
dir.create(paths$out_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------- Parameters -------------------------------------
cluster_use     <- "Cluster_G"  # cluster pour la trajectoire ATAC
n_bins          <- 100
rm_gene_regex   <- c("^LOC", "^NEWGENE", "^RGD")
k_mad_filter    <- 1.5          # seuil initial
anova_qvalue    <- 1e-3
spline_spar     <- 0.45         # valeur par défaut (sera comparée au scan)
scan_spars      <- TRUE         # active le scan spline
k_force_ATAC    <- 4            # K imposé pour ATAC (NA = auto via BIC)
k_force_expr    <- 5
k_force_velo    <- 4
min_seg_len     <- 100
region_id       <- NA           # ex: "chr3-93594270-93594657"
set.seed(42)

# ------------------------------ Helpers ---------------------------------------
minmax_11 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(diff(rng)) || diff(rng) == 0) return(rep(0, length(x)))
  2 * (x - rng[1]) / diff(rng) - 1
}
safe_spline <- function(e, spar) {
  e <- as.numeric(e)
  x <- seq_along(e)
  y <- tryCatch(stats::smooth.spline(x, e, spar = spar)$y, error = function(.) e)
  y
}
filter_by_mad_adaptive <- function(m, k0, k_floor = 0.3, step = 0.2, label = "") {
  if (is.null(m) || length(m) == 0 || nrow(m) == 0) return(m)
  k <- k0
  repeat {
    sdv <- apply(m, 1, sd)
    thr <- median(sdv) + k * mad(sdv)
    sel <- m[sdv > thr, , drop = FALSE]
    if (nrow(sel) > 0 || k <= k_floor) {
      if (nrow(sel) == 0) sel <- m[order(sdv, decreasing = TRUE)[seq_len(min(500, length(sdv)))], , drop = FALSE]
      message(sprintf("[%s] MAD filter: kept %d / %d (k=%.2f).", label, nrow(sel), nrow(m), k))
      return(sel)
    }
    k <- k - step
  }
}

segm_heatmap <- function(mat, obj_label, k_force = NA,
                         heat_cols = viridis::mako(50),
                         cluster_pal = NULL,
                         k_max = 12, min_len = 100) {
  if (is.null(mat) || nrow(mat) == 0) {
    warning(sprintf("%s matrix is empty after filtering.", obj_label))
    return(invisible(NULL))
  }
  plateau <- function(x, thr = .85, min_run = 2) {
    x <- as.numeric(x)
    m <- max(x, na.rm = TRUE)
    if (!is.finite(m) || m == 0) return(which.max(replace(x, !is.finite(x), -Inf)))
    flag <- x / m >= thr
    r <- rle(flag); idx <- which(r$values)
    if (!length(idx)) return(which.max(x))
    st <- cumsum(c(1, head(r$lengths, -1)))[idx]
    en <- st + r$lengths[idx] - 1
    ok <- which((en - st + 1) >= min_run)
    if (!length(ok)) return(which.max(x))
    mid <- (st[ok] + en[ok]) / 2
    as.numeric(mid[which.max(r$lengths[idx][ok])])
  }
  moment  <- apply(mat, 1, plateau)
  ord_idx <- order(moment)
  ord     <- rownames(mat)[ord_idx]
  y       <- as.numeric(moment[ord_idx]); x <- seq_along(ord)
  
  lm0 <- lm(y ~ x)
  bic <- sapply(0:k_max, function(k) {
    mdl <- if (k == 0) lm0 else segmented(lm0, seg.Z = ~x, npsi = k, control = seg.control(display = FALSE))
    BIC(mdl)
  })
  best_k <- which.min(bic) - 1
  if (!is.na(k_force)) best_k <- k_force
  
  if (best_k > 0) {
    mdl_best <- segmented(lm0, seg.Z = ~x, npsi = best_k, control = seg.control(display = FALSE))
    br <- sort(as.numeric(mdl_best$psi[, "Est."]))
    br <- unique(pmin(length(ord) - 1L, pmax(1L, as.integer(round(br)))))
  } else {
    br <- integer(0)
  }
  repeat {
    lens <- diff(c(0, br, length(ord)))
    if (length(lens) > 0 && any(lens < min_seg_len)) {
      drop_at <- which.min(lens)
      if (drop_at == 1) br <- br[-1]
      else if (drop_at == length(lens)) br <- br[-length(br)]
      else br <- br[-(drop_at - 1)]
    } else break
  }
  
  seg_id <- cut(x, c(0, br, length(ord)), labels = FALSE, include.lowest = TRUE)
  K <- length(unique(seg_id))
  pal_fun <- if (is.null(cluster_pal)) function(K) viridis::viridis(K) else cluster_pal
  pal <- pal_fun(K); names(pal) <- as.character(seq_len(K))
  
  hm <- pheatmap(mat[ord, , drop = FALSE],
                 cluster_rows = FALSE, cluster_cols = FALSE,
                 show_rownames = FALSE, show_colnames = FALSE,
                 annotation_row = data.frame(Cluster = factor(seg_id), row.names = ord),
                 annotation_colors = list(Cluster = pal),
                 color = heat_cols,
                 main  = sprintf("%s – %d segments", obj_label, K))
  
  invisible(list(mat_ord = mat[ord, , drop = FALSE],
                 seg_id  = seg_id,
                 breaks  = br,
                 heatmap = hm,
                 K       = K,
                 BIC     = data.frame(k = 0:k_max, BIC = bic)))
}

palette_atac <- function(K) {
  base <- c("#FF9999", "#66B2FF", "#99FF99", "#FFD700", "#FFB266", "#CCCCCC", "#CC99FF")
  if (K <= length(base)) base[1:K] else c(base, viridis(K - length(base)))
}
palette_expr <- function(K) {
  base <- c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#984EA3")
  if (K <= length(base)) base[1:K] else c(base, viridis(K - length(base)))
}
palette_velo <- function(K) {
  base <- c("#EF5350", "#FFB74D", "#FFD54F", "#9CCC65", "#26C6DA", "#42A5F5", "#AB47BC", "#8D6E63")
  if (K <= length(base)) base[1:K] else c(base, viridis(K - length(base)))
}

# =============================================================================
# Load inputs (harmonized) -----------------------------------------------------
# =============================================================================
stopifnot(file.exists(paths$expr_csv), file.exists(paths$vel_csv), file.exists(paths$merged_csv))

# Load Seurat object (prefer existing 'integrated', fall back to RDS or integrated73)
if (!exists("integrated")) {
  if (exists("integrated73")) {
    integrated <- integrated73
  } else {
    stopifnot(file.exists(paths$integrated_rds))
    integrated <- readRDS(paths$integrated_rds)
  }
}
stopifnot(inherits(integrated, "Seurat"))

# Regulatory-time table (single source of truth)
merged_data <- readr::read_csv(paths$merged_csv, col_types = readr::cols())
stopifnot(all(c("Cell", "new_regulatory_time") %in% names(merged_data)))

# Small helpers
minmax_11 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(diff(rng)) || diff(rng) == 0) return(rep(0, length(x)))
  2 * (x - rng[1]) / diff(rng) - 1
}
safe_spline <- function(e, spar) {
  e <- as.numeric(e); x <- seq_along(e)
  y <- tryCatch(stats::smooth.spline(x, e, spar = spar)$y, error = function(.) e)
  y
}
filter_by_mad <- function(m, k) {
  if (is.null(m) || nrow(m) == 0) return(m)
  sdv <- apply(m, 1, sd)
  thr <- median(sdv) + k * mad(sdv)
  m[sdv > thr, , drop = FALSE]
}

# =============================================================================
# ATAC trajectory (from 'integrated') -----------------------------------------
# =============================================================================
message("[3/7] Building ATAC trajectory (Signac::GetTrajectory)…")

# Preconditions
assay_names <- names(integrated@assays)
if (is.null(assay_names) || !"peaks" %in% assay_names) {
  stop("ChromatinAssay 'peaks' not found in `integrated`.")
}

# Subset cluster and map regulatory time (names removed to avoid length/name checks)
objG_chrom <- subset(integrated, subset = cell_Type %in% cluster_use)
stopifnot(ncol(objG_chrom) > 0)
DefaultAssay(objG_chrom) <- "peaks"

rt_lookup <- setNames(merged_data$new_regulatory_time, merged_data$Cell)
cells     <- Cells(objG_chrom)
hit       <- match(cells, names(rt_lookup))

traj_vec <- rep(NA_real_, length(cells))
traj_vec[!is.na(hit)] <- 100 * rt_lookup[hit[!is.na(hit)]]
names(traj_vec) <- NULL
objG_chrom$Trajectory2 <- as.numeric(traj_vec)

# Keep cells with defined trajectory
objG_chrom <- objG_chrom[, !is.na(objG_chrom$Trajectory2)]
stopifnot(ncol(objG_chrom) > 0)

# Compute trajectory
if (!is.function(get0("GetTrajectory"))) {
  stop("GetTrajectory() not found. Load the package that defines it (e.g., Signac).")
}
trajATAC <- suppressMessages(
  GetTrajectory(
    object          = objG_chrom,
    assay           = "peaks",
    trajectory.name = "Trajectory2",
    groupEvery      = 1,
    slot            = "data",
    smoothWindow    = 7,
    log2Norm        = TRUE
  )
)

# Extract matrix
atac_mat <- if ("smoothMat" %in% SummarizedExperiment::assayNames(trajATAC)) {
  SummarizedExperiment::assay(trajATAC, "smoothMat")
} else {
  SummarizedExperiment::assay(trajATAC)
}
if (ncol(atac_mat) != n_bins) {
  message("ATAC matrix has ", ncol(atac_mat), " bins (expected ", n_bins, "). Using actual column count.")
}
built_atac <- TRUE

# =============================================================================
# Velocity (bin by regulatory time) -------------------------------------------
# =============================================================================
message("[4/7] Binning RNA velocity by regulatory time…")
prep_velocity <- function(csv_path, nbins = 100, rm_regex = NULL) {
  stopifnot(file.exists(csv_path))
  vel <- readr::read_csv(csv_path, col_types = readr::cols()) |>
    dplyr::filter(!is.na(new_regulatory_time)) |>
    dplyr::distinct(Cell_ID, .keep_all = TRUE)
  
  meta <- vel |> dplyr::select(Cell_ID, new_regulatory_time)
  expr <- vel |>
    dplyr::select(-dplyr::any_of(c("latent_time", "new_regulatory_time"))) |>
    tibble::column_to_rownames("Cell_ID") |>
    as.matrix() |> t()
  
  meta <- meta[match(colnames(expr), meta$Cell_ID), ]
  bins <- factor(cut(meta$new_regulatory_time, nbins, labels = FALSE, include.lowest = TRUE) - 1,
                 levels = 0:(nbins - 1))
  idx_list <- split(seq_len(ncol(expr)), bins, drop = FALSE)
  vel_bin <- sapply(idx_list, function(idx) {
    if (length(idx) == 0) rep(NA_real_, nrow(expr)) else rowMeans(expr[, idx, drop = FALSE])
  })
  colnames(vel_bin) <- sprintf("T.%d_%d", 0:(nbins - 1), 1:nbins)
  
  if (!is.null(rm_gene_regex)) {
    rm_pat <- paste(rm_gene_regex, collapse = "|")
    vel_bin <- vel_bin[!grepl(rm_pat, rownames(vel_bin)), , drop = FALSE]
  }
  vel_bin
}
vel_mat <- prep_velocity(paths$vel_csv, n_bins, rm_gene_regex)

# =============================================================================
# Expression regtime matrix ----------------------------------------------------
# =============================================================================
message("[5/7] Loading expression regtime matrix…")
expr_mat <- readr::read_csv(paths$expr_csv, col_types = readr::cols()) |>
  dplyr::rename(Gene = 1) |>
  tibble::column_to_rownames("Gene") |>
  as.matrix()
expr_mat <- expr_mat[!grepl(paste(rm_gene_regex, collapse = "|"), rownames(expr_mat)), , drop = FALSE]

# =============================================================================
# Robust variability filter (MAD) ----------------------------------------------
# =============================================================================
message("[6/7] Filtering features by robust variability (MAD)…")
vel_filt  <- filter_by_mad(vel_mat,  k_mad_filter)
expr_filt <- filter_by_mad(expr_mat, k_mad_filter)
atac_filt <- if (isTRUE(built_atac)) filter_by_mad(atac_mat, k_mad_filter) else matrix(numeric(0), 0, 0)

message(sprintf("Kept %d / %d velocity features after MAD filter.",  nrow(vel_filt),  nrow(vel_mat)))
message(sprintf("Kept %d / %d expression genes after MAD filter.", nrow(expr_filt), nrow(expr_mat)))
if (isTRUE(built_atac)) message(sprintf("Kept %d / %d ATAC regions after MAD filter.", nrow(atac_filt), nrow(atac_mat)))

# =============================================================================
# Spline → ANOVA → z-score → min-max ------------------------------------------
# =============================================================================
process_mat <- function(mat, spar = spline_spar, q_cut = anova_qvalue) {
  if (is.null(mat) || nrow(mat) == 0) return(mat)
  
  # 1) smooth each row
  mat_smooth <- t(apply(mat, 1, function(e) safe_spline(e, spar)))
  
  # 2) one-way ANOVA across bins (BH)
  bins <- seq_len(ncol(mat_smooth))
  qv <- apply(mat_smooth, 1, function(e) {
    p <- tryCatch(as.numeric(summary(aov(e ~ bins))[[1]]["bins", "Pr(>F)"]), error = function(err) 1)
    ifelse(is.na(p), 1, p)
  })
  qv <- p.adjust(qv, method = "BH")
  mat_filt <- mat_smooth[qv < q_cut, , drop = FALSE]
  message(nrow(mat_filt), " features passed ANOVA (q < ", q_cut, ")")
  if (nrow(mat_filt) == 0) return(mat_filt)
  
  # 3) row z-score
  mat_scaled <- t(scale(t(mat_filt))); mat_scaled[!is.finite(mat_scaled)] <- 0
  
  # 4) row min-max to [-1,1]
  t(apply(mat_scaled, 1, minmax_11))
}

message("[7/7] Processing matrices (spline → ANOVA → z-score → min-max)…")
velocity_scaled   <- process_mat(vel_filt)
expression_scaled <- process_mat(expr_filt)
chromatin_matrix  <- if (isTRUE(built_atac)) process_mat(atac_filt) else matrix(numeric(0), 0, 0)

# =============================================================================
# Segmentation + heatmaps ------------------------------------------------------
# =============================================================================
res_ATAC <- segm_heatmap(chromatin_matrix,   "ATAC",       k_force = k_force_ATAC,
                         heat_cols = viridis::mako(50),    cluster_pal = palette_atac,
                         min_len = min_seg_len)

res_expr <- segm_heatmap(expression_scaled,  "Expression", k_force = k_force_expr,
                         heat_cols = viridis::plasma(50),  cluster_pal = palette_expr,
                         min_len = min_seg_len)

res_velo <- segm_heatmap(velocity_scaled,    "Velocity",   k_force = k_force_velo,
                         heat_cols = viridis::viridis(50), cluster_pal = palette_velo,
                         min_len = min_seg_len)


# Exports (optionnels)
if (!is.null(res_ATAC))  readr::write_csv(as.data.frame(res_ATAC$BIC),  file.path(paths$out_dir, "BIC_ATAC.csv"))
if (!is.null(res_expr))  readr::write_csv(as.data.frame(res_expr$BIC),  file.path(paths$out_dir, "BIC_Expression.csv"))
if (!is.null(res_velo))  readr::write_csv(as.data.frame(res_velo$BIC),  file.path(paths$out_dir, "BIC_Velocity.csv"))
if (is.matrix(chromatin_matrix) && nrow(chromatin_matrix) > 0)
  write.csv(chromatin_matrix, file.path(paths$out_dir, "ATAC_scaled_matrix.csv"), row.names = TRUE)

message("\n✅ Pipeline complete.")

# ---------------------------------------------------------------------------
# Save segmentation results for downstream ontology analysis (Script 11)
# ---------------------------------------------------------------------------
try({
  save_dir <- paths$out_dir %||% "outputs/heatmaps"
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    list(
      res_expr = res_expr,
      res_ATAC = res_ATAC,
      res_velo = res_velo
    ),
    file = file.path(save_dir, "segmentation_results.rds")
  )
  message("Saved: ", file.path(save_dir, "segmentation_results.rds"))
}, silent = FALSE)
