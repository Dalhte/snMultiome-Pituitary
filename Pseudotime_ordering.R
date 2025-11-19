# ------------------------------------------------------------
# Pseudotime strips (contiguous) — exact AddTrajectory pattern
# Input  : integrated (Seurat), meta: cell_Type, orig.ident
# Needs  : reduction "multimodal_umap", AddTrajectory() available
# Output : Frise_10clusters_nogap.svg (fixed order & colors)
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

integrated <- readRDS('PGintegrated.rds')

# ---- Parameters (keep these consistent with your figure) ----
clusters <- c("G","L","C","S","T","M","FSC","EC","Pe","Le")
reduction_name <- "multimodal_umap"
dims_use <- 1:2
group_by <- "cell_Type"
traj_prefix <- "Trajectory_"
cluster_prefix <- "Cluster_"

out_file <- "Frise_10clusters_nogap.svg"
palette_samples <- c("D2"="#0000FF", "PM"="#00FF00", "PS"="#FF0000", "E"="#FFFF00")

# ---- Minimal sanity checks ----
stopifnot(exists("integrated"), inherits(integrated, "Seurat"))
stopifnot(group_by %in% colnames(integrated@meta.data))
stopifnot("orig.ident" %in% colnames(integrated@meta.data))
stopifnot(reduction_name %in% names(integrated@reductions))
stopifnot(is.function(get0("AddTrajectory")))

# ---- Build rectangles across all clusters (no gaps) ----
all_df_rect <- data.frame()

for (cl in clusters) {
  # 1) Add trajectory directly into integrated (this is the pattern that works for you)
  integrated <- AddTrajectory(
    object     = integrated,
    name       = paste0(traj_prefix, cl),
    trajectory = paste0(cluster_prefix, cl),
    group.by   = group_by,
    reduction  = reduction_name,
    dims       = dims_use,
    use.all    = TRUE
  )
  
  # 2) Subset cells with defined pseudotime for this cluster
  traj_col <- paste0(traj_prefix, cl)
  objTmp <- integrated[, !is.na(integrated[[traj_col]])]
  
  if (ncol(objTmp) == 0) next  # nothing for this cluster
  
  # 3) Ensure orig.ident is present on the subset (keep sample colors)
  objTmp$orig.ident <- integrated$orig.ident
  
  # 4) Build the per-cell table (keep your original extraction idiom)
  df_tmp <- data.frame(
    cell       = colnames(objTmp),
    pseudotime = as.numeric(unlist(objTmp[[traj_col]])),
    orig.ident = objTmp$orig.ident,
    cluster    = cl,
    stringsAsFactors = FALSE
  )
  
  # 5) Sort by pseudotime and compute contiguous x_min/x_max (no holes)
  df_tmp <- df_tmp[order(df_tmp$pseudotime), ]
  df_tmp <- df_tmp %>%
    mutate(
      x_min = (pseudotime + dplyr::lag(pseudotime,  default = pseudotime[1] - 0.001)) / 2,
      x_max = (pseudotime + dplyr::lead(pseudotime, default = pseudotime[length(pseudotime)] + 0.001)) / 2
    )
  
  all_df_rect <- rbind(all_df_rect, df_tmp)
}

stopifnot(nrow(all_df_rect) > 0)

# ---- Fixed order (top to bottom: G … Le) and y mapping ----
all_df_rect$cluster <- factor(all_df_rect$cluster, levels = rev(clusters))
all_df_rect <- all_df_rect %>% mutate(cluster_idx = as.numeric(cluster))

# ---- Plot (same look as your expected figure) ----
p_rect <- ggplot(all_df_rect, aes(fill = orig.ident)) +
  geom_rect(
    aes(
      xmin = x_min, xmax = x_max,
      ymin = cluster_idx - 0.4, ymax = cluster_idx + 0.4
    ),
    color = NA
  ) +
  scale_y_continuous(
    breaks = seq_along(clusters),
    labels = rev(clusters)
  ) +
  scale_fill_manual(values = palette_samples, drop = FALSE) +
  labs(
    title = "Pseudotime strips by cluster (contiguous)",
    x = "Pseudotime",
    y = "Cluster",
    fill = "orig.ident"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(out_file, p_rect, device = "svg", width = 7, height = 5.5, units = "in")
print(p_rect)
