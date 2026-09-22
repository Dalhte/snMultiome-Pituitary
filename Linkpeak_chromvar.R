suppessPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(ggplot2)
})

cat("Loading object...\n")

integrated73 <- readRDS(
"/shared/projects/femcycle/PG_multiome_cluster/results/PGintegrated73.reclustered.rds"
)

cat("Checking peaks assay...\n")

if (!"ATAC" %in% Assays(integrated73)) {
  stop("ATAC assay is missing.")
}

if (!"peaks" %in% Assays(integrated73)) {
  stop("peaks assay is missing.")
}

# peaks was created from ATAC, therefore the per-cell
# ATAC counts/features are also the peaks counts/features.
integrated73$nCount_peaks <- integrated73$nCount_ATAC
integrated73$nFeature_peaks <- integrated73$nFeature_ATAC

stopifnot(
  identical(
    integrated73$nCount_peaks,
    integrated73$nCount_ATAC
  )
)

stopifnot(
  identical(
    integrated73$nFeature_peaks,
    integrated73$nFeature_ATAC
  )
)

cat("nCount_peaks and nFeature_peaks validated.\n")

dir.create(
"/shared/projects/femcycle/PG_multiome_cluster/results/linkpeaks",
showWarnings = FALSE,
recursive = TRUE
)

dir.create(
"/shared/projects/femcycle/PG_multiome_cluster/results/objects",
showWarnings = FALSE,
recursive = TRUE
)

cat("Running LinkPeaks...\n")

integrated73 <- LinkPeaks(
  integrated73,
  peak.assay = "peaks",
  expression.assay = "SoupXRNA",
  peak.slot = "counts",
  expression.slot = "data",
  method = "pearson",
  distance = 5000000,
  min.cells = 10,
  n_sample = 200,
  pvalue_cutoff = 0.05,
  score_cutoff = 0.04,
  verbose = TRUE
)

cat("Extracting links...\n")

links <- Links(integrated73[["peaks"]])

write.csv(
  as.data.frame(links),
  "/shared/projects/femcycle/PG_multiome_cluster/results/linkpeaks/LinkPeaks_results.csv",
  row.names = FALSE
)



meta <- integrated73@meta.data

meta$cluster_pb <- as.character(meta$new_cluster)

meta$cluster_pb[
  meta$cluster_pb %in% c("Cluster_G1", "Cluster_G2")
] <- "Cluster_G"

meta$cluster_pb[
  meta$cluster_pb %in% c("Cluster_L1", "Cluster_L2")
] <- "Cluster_L"

integrated73$cluster_pb <- meta$cluster_pb

integrated73$cluster_pb_short <- sub(
  "^Cluster_",
  "",
  integrated73$cluster_pb
)


desired_cluster_order <- c(
  "S",
  "L",
  "FSC",
  "C",
  "Le",
  "G",
  "T",
  "M",
  "EC",
  "Pe",
  "Pro"
)

integrated73$cluster_pb_short <- factor(
  integrated73$cluster_pb_short,
  levels = desired_cluster_order
)


cluster_color_mapping <- c(
  "S"   = "#F8766D",
  "L"   = "#DB8E00",
  "FSC" = "#AEA200",
  "C"   = "#00BD5C",
  "Le"  = "#00C1A7",
  "G"   = "#00BADE",
  "T"   = "#00A6FF",
  "M"   = "#B385FF",
  "EC"  = "#EF67EB",
  "Pe"  = "#FF63B6",
  "Pro" = "#7F7F7F"
)


levels(integrated73$cluster_pb_short)
table(integrated73$cluster_pb_short)

# ============================================================
# CoveragePlot
# ============================================================

p <- CoveragePlot(
  object = integrated73,
  assay = "peaks",
  region = "Pomc",
  features = "Pomc",
  expression.assay = "SoupXRNA",
  links = TRUE,
  group.by = "cluster_pb_short",
  extend.upstream = 15000,
  extend.downstream = 10000,
  heights = 3
)

p <- p &
  scale_fill_manual(
    values = cluster_color_mapping,
    breaks = desired_cluster_order,
    limits = desired_cluster_order,
    drop = FALSE
  )

p

pdf(
  "/shared/projects/femcycle/PG_multiome_cluster/results/figures/CoveragePlot_Pomc.pdf",
  width = 10,
  height = 8
)

p

dev.off()


p <- CoveragePlot(
  object = integrated73,
  assay = "peaks",
  region = "Nr5a1",
  features = "Nr5a1",
  expression.assay = "SoupXRNA",
  links = TRUE,
  group.by = "cluster_pb_short",
  extend.upstream = 10000,
  extend.downstream = 10000,
  heights = 3
)

p <- p &
  scale_fill_manual(
    values = cluster_color_mapping,
    breaks = desired_cluster_order,
    limits = desired_cluster_order,
    drop = FALSE
  )

p

pdf(
  "/shared/projects/femcycle/PG_multiome_cluster/results/figures/CoveragePlot_Nr5a1.pdf",
  width = 10,
  height = 8
)

p


dev.off()


p <- CoveragePlot(
  object = integrated73,
  assay = "peaks",
  region = "Pou1f1",
  features = "Pou1f1",
  expression.assay = "SoupXRNA",
  links = TRUE,
  group.by = "cluster_pb_short",
  extend.upstream = 10000,
  extend.downstream = 10000,
  heights = 3
)

p <- p &
  scale_fill_manual(
    values = cluster_color_mapping,
    breaks = desired_cluster_order,
    limits = desired_cluster_order,
    drop = FALSE
  )

p

pdf(
  "/shared/projects/femcycle/PG_multiome_cluster/results/figures/CoveragePlot_Pou1f1.pdf",
  width = 10,
  height = 8
)

p

dev.off()

saveRDS(
  integrated73,
  "/shared/projects/femcycle/PG_multiome_cluster/results/objects/PGintegrated73.reclustered.linkpeaks.rds"
)

cat("Finished LinkPeaks\n")