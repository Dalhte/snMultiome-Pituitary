suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(dplyr)
  library(tidyr)
})
cat("Loading object...\n")

integrated73 <- readRDS(
  "/shared/projects/femcycle/PG_multiome_cluster/results/objects/PGintegrated73.reclustered.linkpeaks.rds"
)

cat("Loading filtered LinkPeaks...\n")

linkpeaks <- read.csv(
  "/shared/projects/femcycle/PG_multiome_cluster/results/filtered_linkpeaks/Filtered_LinkPeaks_All.csv",
  stringsAsFactors = FALSE
)

############################################################
# TSS from Signac annotation
############################################################

cat("Building TSS annotation...\n")

ann <- Annotation(integrated73[["peaks"]])

ann_df <- as.data.frame(ann)

tss_df <- ann_df %>%
  mutate(
    strand = as.character(strand)
  ) %>%
  group_by(gene_name) %>%
  summarise(
    TSS_chr = first(as.character(seqnames)),
    strand = first(strand),
    TSS = if(first(strand) == "+"){
      min(start)
    } else {
      max(end)
    },
    .groups = "drop"
  )

cat("Unique genes:", nrow(tss_df), "\n")

############################################################
# Merge LinkPeaks + TSS
############################################################

merged_df <- left_join(
  linkpeaks,
  tss_df,
  by = c("gene" = "gene_name")
)

cat(
  "Links with TSS:",
  sum(!is.na(merged_df$TSS)),
  "/",
  nrow(merged_df),
  "\n"
)

write.csv(
  merged_df,
  "/shared/projects/femcycle/PG_multiome_cluster/results/filtered_linkpeaks/Filtered_LinkPeaks_All_with_TSS.csv",
  row.names = FALSE
)

############################################################
# TAD boundaries
############################################################

cat("Loading TAD boundaries...\n")

tad_boundaries <- read.table(
  "/shared/projects/femcycle/PG_multiome_cluster/data/rat_boundaries_filtered.bed",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

colnames(tad_boundaries) <- c(
  "chr",
  "start",
  "end"
)

tad_boundaries.gr <- makeGRangesFromDataFrame(
  tad_boundaries,
  seqnames.field = "chr",
  start.field = "start",
  end.field = "end",
  ignore.strand = TRUE
)

############################################################
# Recover true peak coordinates
############################################################

peak_coords <- tidyr::extract(
  data.frame(peak = merged_df$peak),
  peak,
  into = c(
    "peak_chr",
    "peak_start",
    "peak_end"
  ),
  regex = "(chr[^-]+)-(\\d+)-(\\d+)",
  convert = TRUE
)

merged_df <- cbind(
  merged_df,
  peak_coords
)

############################################################
# Closest peak edge to TSS
############################################################

merged_df <- merged_df %>%
  mutate(
    dist_start = abs(peak_start - TSS),
    dist_end   = abs(peak_end - TSS),
    candidate_edge = if_else(
      dist_start < dist_end,
      peak_start,
      peak_end
    )
  )
############################################################
# TAD filter function
############################################################

is_edge_separated_by_TAD <- function(
  chr,
  pos1,
  pos2,
  tad_gr
){

  if(is.na(pos1) || is.na(pos2)){
    return(FALSE)
  }

  if(pos1 > pos2){

    tmp <- pos1
    pos1 <- pos2
    pos2 <- tmp

  }

  interval <- GRanges(
    seqnames = chr,
    ranges = IRanges(
      start = pos1,
      end = pos2
    )
  )

  ov <- findOverlaps(
    interval,
    tad_gr
  )

  return(length(ov) > 0)
}

############################################################
# Apply TAD filter
############################################################

cat("Applying TAD filter...\n")

merged_df <- merged_df %>%
  rowwise() %>%
  mutate(
    separated =
      if_else(
  peak_chr == TSS_chr,
  is_edge_separated_by_TAD(
    peak_chr,
    TSS,
    candidate_edge,
    tad_boundaries.gr
  ),
  FALSE
)
  ) %>%
  ungroup()

final_df <- merged_df %>%
  filter(!separated)

cat(
  "Before filter:",
  nrow(merged_df),
  "\n"
)

cat(
  "After filter:",
  nrow(final_df),
  "\n"
)

write.csv(
  final_df,
  "/shared/projects/femcycle/PG_multiome_cluster/results/filtered_linkpeaks/Filtered_LinkPeaks_All_with_TSS_noTAD.csv",
  row.names = FALSE
)

############################################################
# Accessibility heatmap
############################################################

suppressPackageStartupMessages({
  library(Matrix)
  library(matrixStats)
  library(pheatmap)
})

cat("Building accessibility heatmap...\n")

peaks_req <- unique(
  final_df$peak
)

DefaultAssay(integrated73) <- "peaks"

peaks_available <- intersect(
  peaks_req,
  rownames(integrated73[["peaks"]])
)

cat(
  "Peaks found:",
  length(peaks_available),
  "\n"
)

############################################################
# Accessibility matrix
############################################################

access_mat <- GetAssayData(
  integrated73,
  assay = "peaks",
  layer = "data"
)[peaks_available, , drop = FALSE]

############################################################
# Clusters
############################################################

cluster_vec <- integrated73$new_cluster

cluster_vec <- factor(cluster_vec)

design <- Matrix::sparse.model.matrix(
  ~0 + cluster_vec
)

colnames(design) <- levels(cluster_vec)

############################################################
# Cluster averages
############################################################

cells_per_cluster <- Matrix::colSums(design)

avg_by_cluster <- as.matrix(
  access_mat %*% design
)

avg_by_cluster <- sweep(
  avg_by_cluster,
  2,
  cells_per_cluster,
  "/"
)

############################################################
# Z-score
############################################################

row_means <- rowMeans(avg_by_cluster)

row_sds <- matrixStats::rowSds(
  avg_by_cluster
)

z_mat <- sweep(
  avg_by_cluster,
  1,
  row_means,
  "-"
)

z_mat <- sweep(
  z_mat,
  1,
  ifelse(row_sds == 0, 1, row_sds),
  "/"
)

z_mat[
  is.na(z_mat) |
    is.infinite(z_mat)
] <- 0

############################################################
# Save matrix
############################################################

saveRDS(
  z_mat,
  "/shared/projects/femcycle/PG_multiome_cluster/results/objects/DAR_DEG_heatmap_matrix.rds"
)

############################################################
# Heatmap
############################################################

breaks <- seq(
  -4,
  4,
  length.out = 201
)

pal <- colorRampPalette(
  c(
    "black",
    "white",
    "#1aa14e"
  )
)(
  length(breaks) - 1
)

pdf(
  "/shared/projects/femcycle/PG_multiome_cluster/results/figures/Heatmap_DAR_DEG_TADfiltered.pdf",
  width = 10,
  height = 10
)

pheatmap(
  mat = z_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none",
  color = pal,
  breaks = breaks,
  show_rownames = FALSE,
  fontsize_col = 9,
  main = "DAR linked to DEG after TAD filtering",
  clustering_distance_rows = "correlation",
  clustering_method = "ward.D2",
  border_color = NA
)

dev.off()

tiff(
  "/shared/projects/femcycle/PG_multiome_cluster/results/figures/Heatmap_DAR_DEG_TADfiltered.tiff",
  width = 3000,
  height = 3000,
  res = 300
)

pheatmap(
  mat = z_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none",
  color = pal,
  breaks = breaks,
  show_rownames = FALSE,
  fontsize_col = 9,
  main = "DAR linked to DEG after TAD filtering",
  clustering_distance_rows = "correlation",
  clustering_method = "ward.D2",
  border_color = NA
)

dev.off()

cat("Finished\n")