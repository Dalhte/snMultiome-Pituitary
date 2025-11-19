# ================================================================
# Regulatory time (GAM) + regtime exports for expression & velocity
# ================================================================
# Outputs
#   1) outputs/pseudotime_regulatory/merged_regulatory_time.csv
#   2) outputs/regtime/data_regtime.csv                 (optional; needs object RDS + GetTrajectory)
#   3) outputs/regtime/velocity_data_all_regtime.csv    (optional; needs velocity CSV)
# ================================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(mgcv)
  library(patchwork)
  library(scMEGA)
})

# --- Optional packages (used only if available) ---
quiet_require <- function(pkg) suppressWarnings(suppressPackageStartupMessages(require(pkg, character.only = TRUE)))
quiet_require("plotly")
quiet_require("htmlwidgets")
quiet_require("SummarizedExperiment")

# ------------------------------
# 1) Configuration
# ------------------------------
paths <- list(
  pseudotime_csv = "Pseudotime_trajectory_values.csv",
  latenttime_csv = "latenttime_trajectory_values.csv",
  velocity_csv   = "velocity_data_all.csv",   # optional
  object_rds     = "integrated.rds",        # optional (for GetTrajectory)
  out_reg        = "outputs/pseudotime_regulatory",
  out_regtime    = "outputs/regtime"
)

dir.create(paths$out_reg, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$out_regtime, recursive = TRUE, showWarnings = FALSE)

# ------------------------------
# 2) Parsing helpers
# ------------------------------
# Keep your legacy substring rule but add safer fallbacks.
cell_id_from_pseudotime <- function(x, start = 4, end = 15) {
  x <- as.character(x)
  out <- substr(x, start, pmin(end, nchar(x)))
  # if substring collapses to empty for some rows, fall back to full string
  out[nchar(out) == 0] <- x[nchar(out) == 0]
  out
}
cell_type_from_cell <- function(x) sub("^([A-Za-z0-9]+).*", "\\1", as.character(x))

# Color mapping for plots
cell_type_colors <- c(D2 = "#0000FF", PM = "#00FF00", PS = "#FF0000", E = "#FFFF00")

# Deterministic normalizer [0,1]
normalize01 <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

# Deterministic cyclic-GAM fitter:
# - fix basis size (k), scan a fixed grid of sp, single-threaded control
# - choose sp minimizing gcv.ubre (same criterion mgcv prints as “GCV/UBRE”)
fit_gam_cc_deterministic <- function(theta, y, k = 10,
                                     sp_grid = 10^seq(-6, 2, by = 0.25),
                                     score = c("GCV.Cp", "REML")) {
  score <- match.arg(score)
  ctrl  <- mgcv::gam.control(nthreads = 1)
  scans <- lapply(sp_grid, function(spi) {
    m <- mgcv::gam(y ~ s(theta, bs = "cc", k = k),
                   sp = spi, method = score, control = ctrl)
    data.frame(sp = spi, gcv = m$gcv.ubre, stringsAsFactors = FALSE)
  })
  scan_df <- do.call(rbind, scans)
  best_sp <- scan_df$sp[which.min(scan_df$gcv)]
  fit <- mgcv::gam(y ~ s(theta, bs = "cc", k = k),
                   sp = best_sp, method = score, control = ctrl)
  list(fit = fit, best_sp = best_sp, scan = scan_df)
}

# ------------------------------
# 3) Load & harmonize inputs
# ------------------------------
stopifnot(file.exists(paths$pseudotime_csv), file.exists(paths$latenttime_csv))

pt_raw <- readr::read_csv(paths$pseudotime_csv, show_col_types = FALSE)
lt_raw <- readr::read_csv(paths$latenttime_csv, show_col_types = FALSE)

# Detect the cell ID column in each table
detect_cell_col <- function(df, candidates = c("Cell", "cell", "cell_id", "Cell_ID", "barcode")) {
  hit <- intersect(candidates, names(df))
  if (length(hit)) hit[1] else names(df)[1]  # fallback: first col
}
pt_cell_col <- detect_cell_col(pt_raw)
lt_cell_col <- detect_cell_col(lt_raw)

# Detect the pseudotime and latent-time columns
detect_num_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit)) return(hit[1])
  # fallback: first numeric column not used as ID
  num <- names(df)[vapply(df, is.numeric, logical(1))]
  setdiff(num, detect_cell_col(df))[1]
}
pt_val_col <- detect_num_col(pt_raw, c("Trajectory_Value","pseudotime","trajectory_value","PT","pseudotime_value"))
lt_val_col <- detect_num_col(lt_raw, c("latent_time","latentTime","LatentTime","LT","latenttime"))

stopifnot(!is.na(pt_cell_col), !is.na(lt_cell_col), !is.na(pt_val_col), !is.na(lt_val_col))

# Preserve input order with a row id (.pt_row)
pseudotime_data <- pt_raw %>%
  mutate(.pt_row = dplyr::row_number()) %>%
  transmute(
    Cell             = .data[[pt_cell_col]],
    Trajectory_Value = as.numeric(.data[[pt_val_col]]),
    Cell_ID          = cell_id_from_pseudotime(.data[[pt_cell_col]]),
    Cell_Type        = cell_type_from_cell(.data[[pt_cell_col]]),
    .pt_row
  )

latenttime_data <- lt_raw %>%
  transmute(
    Cell_ID     = as.character(.data[[lt_cell_col]]),
    latent_time = as.numeric(.data[[lt_val_col]])
  )

# Keep duplicates if any; align to pseudotime input order
merged_data <- inner_join(pseudotime_data, latenttime_data, by = "Cell_ID") %>%
  arrange(.pt_row)

stopifnot(nrow(merged_data) > 0)

# ------------------------------
# 4) Cyclic GAM → regulatory time (legacy-compatible, deterministic)
# ------------------------------
# Legacy theta: scale by max(Trajectory_Value) only (no min subtraction)
max_tv <- max(merged_data$Trajectory_Value, na.rm = TRUE)
stopifnot(is.finite(max_tv) && max_tv > 0)

merged_data <- merged_data %>%
  mutate(theta = Trajectory_Value * 2 * pi / max_tv)

# Deterministic cyclic GAM on latent_time ~ s(theta, bs="cc")
gam_sel <- fit_gam_cc_deterministic(
  theta   = merged_data$theta,
  y       = merged_data$latent_time,
  k       = 10,
  sp_grid = 10^seq(-6, 2, by = 0.25),
  score   = "GCV.Cp"
)

# Unified time = normalized prediction in [0,1]
uni_raw <- as.numeric(predict(gam_sel$fit, newdata = merged_data))
uni01   <- normalize01(uni_raw)

# Build new_regulatory_time EXACTLY like your legacy:
#  1) arrange by unified_time ONLY
#  2) linear rank: (row_number - 1) / (n - 1)
#  3) restore original file order (.pt_row)
merged_data <- merged_data %>%
  mutate(unified_time = uni01) %>%
  arrange(unified_time) %>%
  mutate(new_regulatory_time = (dplyr::row_number() - 1) / (dplyr::n() - 1)) %>%
  arrange(.pt_row)

# For diagnostics/strips
rank01 <- function(x) {
  x <- as.numeric(x)
  if (length(x) <= 1) return(rep(0, length(x)))
  (rank(x, ties.method = "first") - 1) / (length(x) - 1)
}
merged_data <- merged_data %>%
  mutate(
    new_pseudotime  = rank01(Trajectory_Value),
    new_latent_time = rank01(latent_time)
  )

# ------------------------------
# 5) Diagnostics + strips + friezes (unchanged)
# ------------------------------
# Diagnostic: latent time vs pseudotime with GAM trend
p_scatter <- ggplot(merged_data, aes(Trajectory_Value, latent_time, color = Cell_Type)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = cell_type_colors, drop = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "black") +
  labs(title = "Latent time vs. pseudotime", x = "Pseudotime", y = "Latent time", color = "Group") +
  theme_minimal()
ggsave(file.path(paths$out_reg, "latent_vs_pseudotime.svg"), p_scatter, width = 7, height = 6)

# Stacked strips (PT / LT / RT)
long_df <- merged_data %>%
  select(Cell_Type, new_pseudotime, new_latent_time, new_regulatory_time) %>%
  pivot_longer(cols = -Cell_Type, names_to = "key", values_to = "val") %>%
  mutate(
    strip = recode(key,
                   new_pseudotime      = "PT",
                   new_latent_time     = "LT",
                   new_regulatory_time = "RT"),
    y = recode(strip, PT = 3, LT = 2, RT = 1)
  )

p_strips <- ggplot(long_df, aes(x = val, color = Cell_Type)) +
  geom_segment(aes(xend = val, y = y + 0.4, yend = y - 0.4), linewidth = 1.6) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(breaks = c(1, 2, 3), labels = c("RT", "LT", "PT")) +
  scale_color_manual(values = cell_type_colors, drop = FALSE) +
  labs(title = "Stacked strips of normalized times (PT/LT/RT)", x = "Normalized value [0,1]", y = NULL, color = "Group") +
  theme_minimal() + theme(panel.grid = element_blank())
ggsave(file.path(paths$out_reg, "strips_three_times.svg"), p_strips, width = 7, height = 3)

# Optional: 3D sphere visualization if plotly/htmlwidgets are available
# (display-only; does not affect regulatory time)
if (isTRUE("plotly" %in% .packages(TRUE)) && isTRUE("htmlwidgets" %in% .packages(TRUE))) {
  # Display phi on [0, pi] from latent_time (legacy-style visual)
  phi_disp <- normalize01(merged_data$latent_time) * pi
  theta    <- merged_data$theta
  
  # Fit cyclic GAM: phi_disp ~ s(theta) with proper data frame
  df_fit <- data.frame(theta = theta, phi_disp = phi_disp)
  fit_phi <- mgcv::gam(phi_disp ~ s(theta, bs = "cc", k = 10), data = df_fit)
  
  # Predict smooth curve on a regular theta grid
  theta_seq <- seq(0, 2*pi, length.out = 200)
  gam_df <- tibble::tibble(
    theta = theta_seq,
    phi   = as.numeric(predict(fit_phi, newdata = data.frame(theta = theta_seq)))
  ) %>%
    mutate(
      x = sin(phi) * cos(theta),
      y = sin(phi) * sin(theta),
      z = cos(phi)
    )
  
  # 3D plot
  fig <- plotly::plot_ly() |>
    plotly::add_trace(
      data = merged_data,
      x = ~sin(phi_disp) * cos(theta),
      y = ~sin(phi_disp) * sin(theta),
      z = ~cos(phi_disp),
      type = "scatter3d", mode = "markers",
      marker = list(size = 3.5, opacity = 0.8, color = ~as.factor(Cell_Type)),
      name = "Cells"
    ) |>
    plotly::add_trace(
      data = gam_df, x = ~x, y = ~y, z = ~z,
      type = "scatter3d", mode = "lines",
      line = list(color = "black", width = 8),
      name = "Cyclic GAM"
    ) |>
    plotly::layout(
      title = "Cyclic GAM on the sphere",
      scene = list(xaxis = list(title = "X"), yaxis = list(title = "Y"), zaxis = list(title = "Z"))
    )
  
  htmlwidgets::saveWidget(fig, file = file.path(paths$out_reg, "sphere_cyclic_gam.html"), selfcontained = TRUE)
}


# ------------------------------
# 6) Save merged regulatory table
# ------------------------------
merged_outfile <- file.path(paths$out_reg, "merged_regulatory_time.csv")
readr::write_csv(merged_data, merged_outfile)
message("Saved: ", merged_outfile)

# ================================================================
# 7) OPTIONAL: Export regtime expression matrix (data_regtime.csv)
# ================================================================
# --- Seurat v4 compat layer for Signac::GetTrajectory expecting LayerData (v5 API) ---
`%||%` <- function(a, b) if (!is.null(a)) a else b

if (!exists("LayerData")) {
  LayerData <- function(object, assay = NULL, layer = NULL) {
    # Keep current default unless explicitly set
    if (!is.null(assay)) Seurat::DefaultAssay(object) <- assay
    lyr <- layer %||% "data"
    if (lyr %in% c("data", "scale.data", "counts")) {
      return(Seurat::GetAssayData(object, slot = lyr))
    }
    stop("Unsupported layer: ", lyr,
         " (expected one of: 'counts', 'data', 'scale.data').")
  }
}

export_regtime_expression <- function(object_rds, merged_df, out_csv,
                                      assay_priority = c("SoupXRNA","RNA","SCT"),
                                      traj_name = "Trajectory2",
                                      scale_to_100 = TRUE) {
  stopifnot(file.exists(object_rds))
  obj <- readRDS(object_rds)
  
  # --- colonne temps régul ---
  rt_col <- if ("new_regulatory_time" %in% names(merged_df)) "new_regulatory_time" else
    if ("regulatory_time"      %in% names(merged_df)) "regulatory_time"      else
      stop("No 'new_regulatory_time' or 'regulatory_time' in merged_df.")
  
  # --- clé de mapping : on utilise d'abord 'Cell' (legacy), sinon 'Cell_ID' ---
  if ("Cell" %in% names(merged_df)) {
    key_vec <- as.character(merged_df$Cell)
    val_vec <- as.numeric(merged_df[[rt_col]])
  } else if ("Cell_ID" %in% names(merged_df)) {
    key_vec <- as.character(merged_df$Cell_ID)
    val_vec <- as.numeric(merged_df[[rt_col]])
  } else {
    stop("merged_df needs 'Cell' or 'Cell_ID' to align cells with the Seurat object.")
  }
  
  # --- dictionnaire (première occurrence conserve) ---
  keep_idx <- !duplicated(key_vec)
  map_keys <- key_vec[keep_idx]
  map_vals <- val_vec[keep_idx]
  names(map_vals) <- map_keys
  
  # --- appariement sur colnames(obj) ---
  cells_obj <- colnames(obj)
  idx <- match(cells_obj, names(map_vals))
  has_match <- !is.na(idx)
  n_match <- sum(has_match)
  if (n_match == 0) stop("No overlap between Seurat colnames(obj) and merged_df keys.")
  
  # --- sous-échantillonne l'objet AVANT d'écrire la méta ---
  obj <- obj[, has_match, drop = FALSE]
  idx2 <- idx[has_match]
  traj2 <- map_vals[idx2]
  if (isTRUE(scale_to_100)) traj2 <- 100 * traj2
  # IMPORTANT: vecteur non nommé, longueur == ncol(obj)
  names(traj2) <- NULL
  obj[[traj_name]] <- traj2
  
  # --- choisir l'assay RNA disponible ---
  assays_available <- names(obj@assays)
  assay_use <- intersect(assay_priority, assays_available)
  if (!length(assay_use)) {
    stop("No RNA assay found. Available assays: ", paste(assays_available, collapse = ", "))
  }
  assay_use <- assay_use[1]
  Seurat::DefaultAssay(obj) <- assay_use
  
  # --- GetTrajectory ---
  if (!requireNamespace("Signac", quietly = TRUE)) stop("Package 'Signac' required for GetTrajectory().")
  if (!is.function(get0("GetTrajectory"))) stop("GetTrajectory() not found in the session.")
  traj <- suppressMessages(
    GetTrajectory(
      object          = obj,
      assay           = assay_use,
      trajectory.name = traj_name,
      groupEvery      = 1,
      slot            = "data",
      smoothWindow    = 7,
      log2Norm        = TRUE
    )
  )
  
  # --- extraire la matrice ---
  get_mat <- function(se) {
    if (inherits(se, "SummarizedExperiment")) {
      anns <- SummarizedExperiment::assayNames(se)
      if ("smoothMat" %in% anns) SummarizedExperiment::assay(se, "smoothMat") else SummarizedExperiment::assay(se)
    } else {
      as.matrix(se)
    }
  }
  mat <- get_mat(traj)
  
  # --- filtre gènes et export ---
  row_keep <- !grepl("^(LOC|NEWGENE|RGD)", rownames(mat))
  mat <- mat[row_keep, , drop = FALSE]
  out_tbl <- tibble::tibble(gene = rownames(mat)) |>
    dplyr::bind_cols(as.data.frame(mat, check.names = FALSE))
  readr::write_csv(out_tbl, out_csv)
  message("Saved: ", out_csv, " (", nrow(mat), " genes x ", ncol(mat), " bins; assay=", assay_use, ")")
  invisible(TRUE)
}

# --- appel identique ---
regtime_csv <- file.path(paths$out_regtime, "data_regtime.csv")
if (file.exists(paths$object_rds)) {
  suppressPackageStartupMessages(library(Signac))
  try(export_regtime_expression(
    object_rds   = paths$object_rds,
    merged_df    = merged_data,
    out_csv      = regtime_csv,
    assay_priority = c("SoupXRNA","RNA","SCT"),
    traj_name    = "Trajectory2",
    scale_to_100 = TRUE
  ), silent = FALSE)
} else {
  message("Skip data_regtime.csv (object RDS not found at ", paths$object_rds, ").")
}

# ================================================================
# 8) OPTIONAL: Merge regulatory time into velocity table
# ================================================================
regtime_vel_csv <- file.path(paths$out_regtime, "velocity_data_all_regtime.csv")
if (file.exists(paths$velocity_csv)) {
  velocity_raw <- readr::read_csv(paths$velocity_csv, show_col_types = FALSE)
  if (!"Cell_ID" %in% names(velocity_raw)) names(velocity_raw)[1] <- "Cell_ID"
  
  if (!"Cell_ID" %in% names(merged_data)) {
    merged_data <- merged_data %>%
      mutate(Cell_ID = if ("Cell_ID" %in% names(.)) Cell_ID else cell_id_from_pseudotime(Cell))
  }
  
  vel_joined <- velocity_raw %>%
    left_join(merged_data %>% select(Cell_ID, new_regulatory_time), by = "Cell_ID") %>%
    filter(!is.na(new_regulatory_time))
  
  readr::write_csv(vel_joined, regtime_vel_csv)
  message("Saved: ", regtime_vel_csv, " (", nrow(vel_joined), " cells)")
} else {
  message("Skip velocity_data_all_regtime.csv (velocity CSV not found at ", paths$velocity_csv, ").")
}

# ================================================================
# Friezes ordered by their own axis (PT / LT / RT) — one tile per cell
# ================================================================

stopifnot(all(c("new_pseudotime","new_latent_time","new_regulatory_time","Cell_Type") %in% names(merged_data)))

library(scales)

# Long format: one row per cell × time scale
fr_long <- merged_data %>%
  transmute(
    Cell_Type,
    PT = new_pseudotime,
    LT = new_latent_time,
    RT = new_regulatory_time
  ) %>%
  pivot_longer(PT:RT, names_to = "strip", values_to = "val")

# Rank within strip -> per-cell x in [0,1] and per-strip width
fr_long <- fr_long %>%
  group_by(strip) %>%
  arrange(val, .by_group = TRUE) %>%
  mutate(
    n_strip = n(),
    rank    = row_number(),
    x       = (rank - 0.5) / n_strip,  # center of each tile
    w       = 1 / n_strip              # tile width
  ) %>%
  ungroup() %>%
  mutate(
    strip_lab = factor(strip,
                       levels = c("PT","LT","RT"),
                       labels = c("Pseudotime","Latent time","Regulatory time")
    )
  )

# Fixed palette (same as earlier scripts)
pal_stage <- c(D2 = "#0000FF", PM = "#00FF00", PS = "#FF0000", E = "#FFFF00")

p_friezes <- ggplot(fr_long, aes(x = x, y = strip_lab, fill = Cell_Type)) +
  geom_tile(aes(width = w, height = 0.9), color = NA) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, .25, .5, .75, 1),
    labels = label_percent(accuracy = 1)
  ) +
  scale_fill_manual(values = pal_stage, drop = FALSE, name = "Stage") +
  labs(
    title = "Friezes by own ordering (PT / LT / RT)",
    x = "0 … 100%",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 11)
  )

print(p_friezes)
ggsave(file.path(paths$out_reg, "friezes_by_own_order.svg"),
       p_friezes, width = 8, height = 3, units = "in")
