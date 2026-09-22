suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(rtracklayer)
  library(Matrix)
})

set.seed(1234)

cat("Loading object...\n")

obj <- readRDS(
  "/shared/projects/femcycle/PG_multiome_cluster/data/PGintegrated73.rebuilt-harm_cluster2new_removed.rds"
)


################################################
# CREATE PEAKS ASSAY FROM ATAC
################################################

cat("\n")
cat("================================================\n")
cat("Creating peaks assay from ATAC...\n")
cat("================================================\n")

if (!"ATAC" %in% Assays(obj)) {
  stop("ERROR: ATAC assay is missing from object.")
}

if ("peaks" %in% Assays(obj)) {
  stop(
    "ERROR: peaks assay already exists. ",
    "Refusing to overwrite it."
  )
}

obj[["peaks"]] <- obj[["ATAC"]]

cat(
  "ATAC features:",
  nrow(obj[["ATAC"]]),
  "\n"
)

cat(
  "peaks features:",
  nrow(obj[["peaks"]]),
  "\n"
)

stopifnot(
  identical(
    rownames(obj[["ATAC"]]),
    rownames(obj[["peaks"]])
  )
)

stopifnot(
  identical(
    colnames(obj[["ATAC"]]),
    colnames(obj[["peaks"]])
  )
)

cat("ATAC -> peaks copy validated.\n")



################################################
# FIX PEAK ANNOTATION
################################################

cat("\n")
cat("================================================\n")
cat("Fixing peaks annotation...\n")
cat("================================================\n")

DefaultAssay(obj) <- "peaks"

gtf.path <- "/shared/projects/femcycle/PG_multiome_cluster/data/genes/genes.gtf"

if (!file.exists(gtf.path)) {
  stop(
    "ERROR: Cell Ranger GTF not found:\n",
    gtf.path
  )
}

cat("Importing Cell Ranger GTF...\n")

annotation.correct <- rtracklayer::import(
  gtf.path
)

cat(
  "Imported GTF:",
  length(annotation.correct),
  "ranges\n"
)

################################################
# REQUIRED SIGNAC METADATA
################################################

cat("Preparing Signac annotation metadata...\n")

# Cell Ranger GTF contains "original_biotype"
# Signac requires "gene_biotype"

if (
  !"gene_biotype" %in% colnames(mcols(annotation.correct))
) {

  if (
    "original_biotype" %in%
    colnames(mcols(annotation.correct))
  ) {

    annotation.correct$gene_biotype <-
      annotation.correct$original_biotype

  } else {

    stop(
      "ERROR: neither gene_biotype nor original_biotype ",
      "exists in the GTF."
    )
  }
}

################################################
# Check required columns
################################################

required.columns <- c(
  "gene_name",
  "gene_id",
  "gene_biotype",
  "type"
)

missing.columns <- setdiff(
  required.columns,
  colnames(mcols(annotation.correct))
)

if (length(missing.columns) > 0) {

  stop(
    "ERROR: missing required annotation columns: ",
    paste(
      missing.columns,
      collapse = ", "
    )
  )
}

cat(
  "Required columns present:",
  paste(
    required.columns,
    collapse = ", "
  ),
  "\n"
)

################################################
# Keep relevant annotation records
################################################

annotation.correct <- annotation.correct[
  annotation.correct$type %in%
    c(
      "gene",
      "transcript",
      "exon"
    )
]

################################################
# Set genome
################################################

genome(annotation.correct) <- "mRatBN7.2"

################################################
# Replace corrupted annotation
################################################

cat("Replacing corrupted peaks annotation...\n")

Annotation(
  obj[["peaks"]]
) <- annotation.correct

################################################
# VALIDATE
################################################

cat("\n")
cat("================================================\n")
cat("Validating corrected annotation...\n")
cat("================================================\n")

ann.check <- Annotation(
  obj[["peaks"]]
)

genes.check <- c(
  "Tyrobp",
  "Gh1",
  "Pax7",
  "Lhb",
  "Mki67"
)

expected.chr <- c(
  "chr1",
  "chr10",
  "chr5",
  "chr1",
  "chr1"
)

for (i in seq_along(genes.check)) {

  g <- genes.check[i]
  expected <- expected.chr[i]

  x <- ann.check[
    !is.na(ann.check$gene_name) &
    ann.check$gene_name == g
  ]

  if (length(x) == 0) {

    stop(
      "ERROR: gene missing from corrected annotation: ",
      g
    )
  }

  observed <- unique(
    as.character(
      seqnames(x)
    )
  )

  cat(
    g,
    ": expected",
    expected,
    "| observed:",
    paste(
      observed,
      collapse = ", "
    ),
    "\n"
  )

  if (!(expected %in% observed)) {

    stop(
      "ERROR: incorrect chromosome for ",
      g,
      ". Expected ",
      expected,
      " but found ",
      paste(
        observed,
        collapse = ", "
      )
    )
  }
}

cat("\nAnnotation validation PASSED.\n")

################################################
# GeneActivity
################################################

cat("\n")
cat("================================================\n")
cat("Building GeneActivity...\n")
cat("================================================\n")

DefaultAssay(obj) <- "peaks"

################################################
# Check annotation
################################################

ann.check <- Annotation(
  obj[["peaks"]]
)

cat(
  "Annotation:",
  length(ann.check),
  "ranges\n"
)

cat(
  "Genes:",
  length(
    unique(
      na.omit(
        ann.check$gene_name
      )
    )
  ),
  "\n"
)

cat("\nGene biotypes present:\n")

print(
  sort(
    table(
      ann.check$gene_biotype,
      useNA = "always"
    ),
    decreasing = TRUE
  )
)

################################################
# Validate required genes
################################################

genes.check <- c(
  "Tyrobp",
  "Gh1",
  "Pax7",
  "Lhb",
  "Mki67"
)

expected.chr <- c(
  "chr1",
  "chr10",
  "chr5",
  "chr1",
  "chr1"
)

for (i in seq_along(genes.check)) {

  g <- genes.check[i]
  expected <- expected.chr[i]

  x <- ann.check[
    !is.na(ann.check$gene_name) &
    ann.check$gene_name == g
  ]

  if (length(x) == 0) {
    stop(
      "ERROR: gene missing from annotation: ",
      g
    )
  }

  observed <- unique(
    as.character(
      seqnames(x)
    )
  )

  cat(
    g,
    ": expected",
    expected,
    "| observed:",
    paste(
      observed,
      collapse = ", "
    ),
    "\n"
  )

  if (!(expected %in% observed)) {

    stop(
      "ERROR: incorrect chromosome for ",
      g,
      ". Expected ",
      expected,
      " but found ",
      paste(
        observed,
        collapse = ", "
      )
    )
  }
}

cat("\nAnnotation validation PASSED.\n")

################################################
# GeneActivity
################################################

cat("\n")
cat("Running GeneActivity with all annotation biotypes...\n")

gene.activity <- GeneActivity(
  object = obj,
  biotypes = NULL,
  process_n = 10000,
  verbose = TRUE
)

################################################
# Validate GeneActivity
################################################

cat("\n")
cat("================================================\n")
cat("GeneActivity validation\n")
cat("================================================\n")

for (g in genes.check) {

  if (!(g %in% rownames(gene.activity))) {

    stop(
      "ERROR: ",
      g,
      " is absent from GeneActivity matrix."
    )
  }

  v <- gene.activity[g, ]

  nz <- sum(v > 0)
  total <- sum(v)

  cat(
    g,
    ":",
    nz,
    "cells > 0 | total activity =",
    total,
    "\n"
  )

  if (nz == 0) {

    stop(
      "ERROR: GeneActivity for ",
      g,
      " is still zero in every cell."
    )
  }
}

cat("\nGeneActivity validation PASSED.\n")

################################################
# Create GeneActivity assay
################################################

obj[["GeneActivity"]] <- CreateAssayObject(
  counts = gene.activity
)

obj <- NormalizeData(
  object = obj,
  assay = "GeneActivity",
  normalization.method = "LogNormalize",
  scale.factor = median(
    Matrix::colSums(
      GetAssayData(
        obj,
        assay = "GeneActivity",
        layer = "counts"
      )
    )
  )
)

DefaultAssay(obj) <- "SoupXRNA"

cat(
  "\nGeneActivity:",
  nrow(obj[["GeneActivity"]]),
  "genes x",
  ncol(obj[["GeneActivity"]]),
  "cells\n"
)

gc()

cat("Cells:", ncol(obj), "\n")



################################################
# Rebuild WNN
################################################

obj <- FindMultiModalNeighbors(
  object = obj,
  reduction.list = list("harmony_pca", "harmony_lsi"),
  dims.list = list(2:30, 2:30),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

################################################
# UMAP
################################################

obj <- RunUMAP(
  object = obj,
  nn.name = "weighted.nn",
  reduction.name = "multimodal_umap_2D",
  assay = "RNA",
  verbose = TRUE,
  return.neighbor = TRUE
)

################################################
# Clustering
################################################

obj <- FindClusters(
  obj,
  resolution = 0.45,
  graph = "wknn",
  verbose = TRUE
)

################################################
# New cluster annotation
################################################

meta <- obj@meta.data

meta$new_cluster <- as.character(meta$cluster_id)

meta$new_cluster[
  meta$seurat_clusters %in% c(6, 0)
] <- "Cluster_L1"

meta$new_cluster[
  meta$seurat_clusters %in% c(2)
] <- "Cluster_L2"

meta$new_cluster[
  meta$seurat_clusters %in% c(5)
] <- "Cluster_G1"

meta$new_cluster[
  meta$seurat_clusters %in% c(10)
] <- "Cluster_G2"

meta$new_cluster[
  meta$seurat_clusters %in% c(1)
] <- "Cluster_S"

meta$new_cluster[
  meta$seurat_clusters %in% c(9)
] <- "Cluster_T"

meta$new_cluster[
  meta$seurat_clusters %in% c(3)
] <- "Cluster_FSC"

meta$new_cluster[
  meta$seurat_clusters %in% c(12)
] <- "Cluster_Pe"

meta$new_cluster[
  meta$seurat_clusters %in% c(7)
] <- "Cluster_EC"

meta$new_cluster[
  meta$seurat_clusters %in% c(11)
] <- "Cluster_M"

meta$new_cluster[
  meta$seurat_clusters %in% c(4)
] <- "Cluster_C"

meta$new_cluster[
  meta$seurat_clusters %in% c(8)
] <- "Cluster_Le"

meta$new_cluster[
  meta$seurat_clusters %in% c(13)
] <- "Cluster_Pro"

obj$new_cluster <- meta$new_cluster

################################################
# UMAP figures
################################################

cat("Exporting UMAP figures...\n")

DefaultAssay(obj) <- "SoupXRNA"

p1 <- DimPlot(
  obj,
  reduction = "multimodal_umap_2D",
  group.by = "new_cluster",
  label = FALSE
)

p2 <- DimPlot(
  obj,
  reduction = "multimodal_umap_2D",
  group.by = "new_cluster",
  label = TRUE,
  repel = TRUE
)

pdf(
  "/shared/projects/femcycle/PG_multiome_cluster/results/figures/UMAP_clusters.pdf",
  width = 8,
  height = 7
)

print(p1)
print(p2)

dev.off()

tiff(
  "/shared/projects/femcycle/PG_multiome_cluster/results/figures/UMAP_clusters.tiff",
  width = 3000,
  height = 2500,
  res = 300
)

print(p2)

dev.off()

################################################
# Publication figures
################################################

dir.create(
  "/shared/projects/femcycle/PG_multiome_cluster/results/figures",
  showWarnings = FALSE,
  recursive = TRUE
)

cluster_colors <- c(
  "Cluster_S"   = "#F8766D",
  "Cluster_L1"  = "#DB8E00",
  "Cluster_L2"  = "#DB8E00",
  "Cluster_FSC" = "#AEA200",
  "Cluster_C"   = "#00BD5C",
  "Cluster_Le"  = "#00C1A7",
  "Cluster_G1"  = "#00BADE",
  "Cluster_G2"  = "#00BADE",
  "Cluster_T"   = "#00A6FF",
  "Cluster_M"   = "#B385FF",
  "Cluster_EC"  = "#EF67EB",
  "Cluster_Pe"  = "#FF63B6",
  "Cluster_Pro" = "#7F7F7F"
)

################################################
# UMAP publication
################################################

p_umap <- DimPlot(
  obj,
  reduction = "multimodal_umap_2D",
  group.by = "new_cluster",
  label = FALSE,
  pt.size = 0.05,
  cols = cluster_colors
)

pdf(
  "/shared/projects/femcycle/PG_multiome_cluster/results/figures/UMAP_clusters.pdf",
  width = 8,
  height = 6
)

print(p_umap)

dev.off()

tiff(
  "/shared/projects/femcycle/PG_multiome_cluster/results/figures/UMAP_clusters.tiff",
  width = 3000,
  height = 2500,
  res = 300
)

print(p_umap)

dev.off()

################################################
# FeaturePlots SoupXRNA
################################################

genes_of_interest <- c(
  "Emcn","Plvap","Pde5a","Pdgfrb","Tyrobp",
  "Arhgap15","Adamts9","Rfx4","S100b",
  "Prl","S100g","Pou1f1","Gh1","Tshb",
  "Cga","Pomc","Neurod1","Pax7",
  "Lhb","Fshb","Gnrhr","Nr5a1", "Mki67"
)

DefaultAssay(obj) <- "SoupXRNA"

p_feat <- FeaturePlot(
  obj,
  reduction = "multimodal_umap_2D",
  features = genes_of_interest,
  max.cutoff = "q90",
  min.cutoff = 0,
  pt.size = 0.05,
  order = TRUE
)

pdf(
  "/shared/projects/femcycle/PG_multiome_cluster/results/figures/FeaturePlot_SoupXRNA.pdf",
  width = 16,
  height = 12
)

print(p_feat)

dev.off()

tiff(
  "/shared/projects/femcycle/PG_multiome_cluster/results/figures/FeaturePlot_SoupXRNA.tiff",
  width = 6000,
  height = 4500,
  res = 300
)

print(p_feat)

dev.off()

################################################
# FeaturePlots GeneActivity
################################################

DefaultAssay(obj) <- "GeneActivity"

p_feat <- FeaturePlot(
  obj,
  reduction = "multimodal_umap_2D",
  features = genes_of_interest,
  max.cutoff = "q90",
  min.cutoff = 0,
  pt.size = 0.05,
  order = TRUE
)

pdf(
  "/shared/projects/femcycle/PG_multiome_cluster/results/figures/FeaturePlot_GeneActivity.pdf",
  width = 16,
  height = 12
)

print(p_feat)

dev.off()

tiff(
  "/shared/projects/femcycle/PG_multiome_cluster/results/figures/FeaturePlot_GeneActivity.tiff",
  width = 6000,
  height = 4500,
  res = 300
)

print(p_feat)

dev.off()

DefaultAssay(obj) <- "SoupXRNA"


################################################
# Save object
################################################

saveRDS(
  obj,
  file = "/shared/projects/femcycle/PG_multiome_cluster/results/PGintegrated73.reclustered.rds"
)

write.csv(
  table(obj$new_cluster),
  "/shared/projects/femcycle/PG_multiome_cluster/results/cluster_counts.csv"
)

cat("Done\n")