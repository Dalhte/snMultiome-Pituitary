###################################################################################################################################################
#                                                                                                                                                 #
#                                                      Multimodal integration of the 4 samples                                                    #
#                                                                                                                                                 #
###################################################################################################################################################


#Load libraries

.libPaths('your/path/to/script')

library(Seurat)
library(Signac)
options(Seurat.object.assay.version = 'v5')
library(GenomicRanges)
library(ggplot2)
library(future)
library(reticulate)
library(umap)
library(tidyverse)
library(SoupX)
library(RColorBrewer)
library(nucleR)
library(future)

library(EnsDb.Rnorvegicus.v110)
library(BSgenome.Rnorvegicus.NCBI.rn7.2)
library(universalmotif)
library(TFBSTools)


set.seed(1234)
setwd('/your/path/to/directory')
plan("multisession", workers = 10)
options(future.globals.maxSize = 50000 * 1024^2) 


###################################################################################################################################################
#                                                      Ambiant RNA removal with soupX                                                             #
###################################################################################################################################################


#PG2
# Load data
filt.matrixPG2 <- Read10X_h5("path/to/PG2/outs/filtered_feature_bc_matrix.h5",use.names = T)
raw.matrixPG2  <- Read10X_h5("path/to/PG2/outs/raw_feature_bc_matrix.h5",use.names = T)
cts.raw.matrixPG2 <- raw.matrixPG2$`Gene Expression`
cts.filt.matrixPG2 <- filt.matrixPG2$`Gene Expression`

#Soupx
#quick clustering
sratPG2  <- CreateSeuratObject(counts = cts.filt.matrixPG2)
Soup.channelPG2  <- SoupChannel(cts.raw.matrixPG2, cts.filt.matrixPG2)
sratPG2    <- SCTransform(sratPG2, verbose = FALSE)
sratPG2    <- RunPCA(sratPG2, verbose = F)
sratPG2    <- RunUMAP(sratPG2, dims = 1:15, verbose = F)
sratPG2    <- FindNeighbors(sratPG2, dims = 1:15, verbose = F)
sratPG2    <- FindClusters(sratPG2, verbose = T)

# add clustering to the channel using setClusters. setDR is useful for visualizations.
metaPG2    <- sratPG2@meta.data
umapPG2    <- sratPG2@reductions$umap@cell.embeddings
Soup.channelPG2  <- setClusters(Soup.channelPG2, metaPG2$seurat_clusters)
Soup.channelPG2  <- setDR(Soup.channelPG2, umapPG2)

#run the main SoupX function, calculating ambient RNA profile.
Soup.channelPG2  <- autoEstCont(Soup.channelPG2)
Soup.channelPG2  <- setContaminationFraction(Soup.channelPG2, 0.20, forceAccept = FALSE)

#roundToInt option to make sure we output integer matrix.
adj.matrixPG2  <- adjustCounts(Soup.channelPG2, roundToInt = T)

# write the directory with corrected read counts
DropletUtils:::write10xCounts("path/to/SoupX/SoupX_PG2_RN73", adj.matrixPG2)
saveRDS(adj.matrixPG2,"path/to/SoupX/soupx_corrected_countsPG2RN73.rds")

#PG24
# Load data
filt.matrixPG24 <- Read10X_h5("path/to/PG24/outs/filtered_feature_bc_matrix.h5",use.names = T)
raw.matrixPG24  <- Read10X_h5("path/to/PG24/outs/raw_feature_bc_matrix.h5",use.names = T)
cts.raw.matrixPG24 <- raw.matrixPG24$`Gene Expression`
cts.filt.matrixPG24 <- filt.matrixPG24$`Gene Expression`

#Soupx
#quick clustering
sratPG24  <- CreateSeuratObject(counts = cts.filt.matrixPG24)
Soup.channelPG24  <- SoupChannel(cts.raw.matrixPG24, cts.filt.matrixPG24)
sratPG24    <- SCTransform(sratPG24, verbose = FALSE)
sratPG24    <- RunPCA(sratPG24, verbose = F)
sratPG24    <- RunUMAP(sratPG24, dims = 1:15, verbose = F)
sratPG24    <- FindNeighbors(sratPG24, dims = 1:15, verbose = F)
sratPG24    <- FindClusters(sratPG24, verbose = T)

# add clustering to the channel using setClusters. setDR is useful for visualizations.
metaPG24    <- sratPG24@meta.data
umapPG24    <- sratPG24@reductions$umap@cell.embeddings
Soup.channelPG24  <- setClusters(Soup.channelPG24, metaPG24$seurat_clusters)
Soup.channelPG24  <- setDR(Soup.channelPG24, umapPG24)

#run the main SoupX function, calculating ambient RNA profile.
Soup.channelPG24  <- autoEstCont(Soup.channelPG24)
Soup.channelPG24  <- setContaminationFraction(Soup.channelPG24, 0.20, forceAccept = FALSE)

#roundToInt option to make sure we output integer matrix.
adj.matrixPG24  <- adjustCounts(Soup.channelPG24, roundToInt = T)

# write the directory with corrected read counts.
DropletUtils:::write10xCounts("path/to/SoupX/soupX_PG24_RN73", adj.matrixPG24)
saveRDS(adj.matrixPG24,"path/to/SoupX/soupx_corrected_countsPG24RN73.rds")

#PG6
# Load data 
filt.matrixPG6 <- Read10X_h5("path/to/PG6/outs/filtered_feature_bc_matrix.h5",use.names = T)
raw.matrixPG6  <- Read10X_h5("path/to/PG6/outs/raw_feature_bc_matrix.h5",use.names = T)
cts.raw.matrixPG6 <- raw.matrixPG6$`Gene Expression`
cts.filt.matrixPG6 <- filt.matrixPG6$`Gene Expression`

#Soupx
#quick clustering
sratPG6  <- CreateSeuratObject(counts = cts.filt.matrixPG6)
Soup.channelPG6  <- SoupChannel(cts.raw.matrixPG6, cts.filt.matrixPG6)
sratPG6    <- SCTransform(sratPG6, verbose = FALSE)
sratPG6    <- RunPCA(sratPG6, verbose = F)
sratPG6    <- RunUMAP(sratPG6, dims = 1:15, verbose = F)
sratPG6    <- FindNeighbors(sratPG6, dims = 1:15, verbose = F)
sratPG6    <- FindClusters(sratPG6, verbose = T)

# add clustering to the channel using setClusters. setDR is useful for visualizations.
metaPG6    <- sratPG6@meta.data
umapPG6    <- sratPG6@reductions$umap@cell.embeddings
Soup.channelPG6  <- setClusters(Soup.channelPG6, metaPG6$seurat_clusters)
Soup.channelPG6  <- setDR(Soup.channelPG6, umapPG6)

#run the main SoupX function, calculating ambient RNA profile.
Soup.channelPG6  <- autoEstCont(Soup.channelPG6)
Soup.channelPG6  <- setContaminationFraction(Soup.channelPG6, 0.10, forceAccept = FALSE)

#roundToInt option to make sure we output integer matrix.
adj.matrixPG6  <- adjustCounts(Soup.channelPG6, roundToInt = T)

# write the directory with corrected read counts.
DropletUtils:::write10xCounts("path/to/SoupX/soupX_PG6_RN73", adj.matrixPG6)
saveRDS(adj.matrixPG6,"path/to/SoupX/soupx_corrected_countsPG6RN73.rds")

#PG13
#Load data
filt.matrixPG13 <- Read10X_h5("path/to/PG13/outs/filtered_feature_bc_matrix.h5",use.names = T)
raw.matrixPG13  <- Read10X_h5("path/to/PG13/outs/raw_feature_bc_matrix.h5",use.names = T)
cts.raw.matrixPG13 <- raw.matrixPG13$`Gene Expression`
cts.filt.matrixPG13 <- filt.matrixPG13$`Gene Expression`


#Soupx
#quick clustering
sratPG13  <- CreateSeuratObject(counts = cts.filt.matrixPG13)
Soup.channelPG13  <- SoupChannel(cts.raw.matrixPG13, cts.filt.matrixPG13)
sratPG13    <- SCTransform(sratPG13, verbose = FALSE)
sratPG13    <- RunPCA(sratPG13, verbose = F)
sratPG13    <- RunUMAP(sratPG13, dims = 1:15, verbose = F)
sratPG13    <- FindNeighbors(sratPG13, dims = 1:15, verbose = F)
sratPG13    <- FindClusters(sratPG13, verbose = T)

# add clustering to the channel using setClusters. setDR is useful for visualizations.
metaPG13    <- sratPG13@meta.data
umapPG13    <- sratPG13@reductions$umap@cell.embeddings
Soup.channelPG13  <- setClusters(Soup.channelPG13, metaPG13$seurat_clusters)
Soup.channelPG13  <- setDR(Soup.channelPG13, umapPG13)

#run the main SoupX function, calculating ambient RNA profile.
Soup.channelPG13  <- autoEstCont(Soup.channelPG13)
Soup.channelPG13  <- setContaminationFraction(Soup.channelPG13, 0.12, forceAccept = FALSE)

#roundToInt option to make sure we output integer matrix.
adj.matrixPG13  <- adjustCounts(Soup.channelPG13, roundToInt = T)

# write the directory with corrected read counts.
DropletUtils:::write10xCounts("path/to/SoupX/soupX_PG13_RN73", adj.matrixPG13)
saveRDS(adj.matrixPG13,"path/to/SoupX/soupx_corrected_countsPG13RN73.rds")


###################################################################################################################################################
#                                                      Preparation of chromatine data                                                             #
###################################################################################################################################################



### get gene annotations for rn7
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Rnorvegicus.v110)
seqlevelsStyle(annotation) <- "UCSC"



### Creation of a shared fragment database

peaks.PG2 <- read.table(
  file = "PG2RN73/outs/atac_peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.PG6 <- read.table(
  file = "PG6RN73/outs/atac_peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.PG24 <- read.table(
  file = "PG24RN73/outs/atac_peaks.bed",
  col.names = c("chr", "start", "end")
)
peaks.PG13 <- read.table(
  file = "PG13RN73/outs/atac_peaks.bed",
  col.names = c("chr", "start", "end")
)

# convert to genomic ranges
gr.PG2 <- makeGRangesFromDataFrame(peaks.PG2)
gr.PG6 <- makeGRangesFromDataFrame(peaks.PG6)
gr.PG24 <- makeGRangesFromDataFrame(peaks.PG24)
gr.PG13 <- makeGRangesFromDataFrame(peaks.PG13)


# Create a unified set of peaks to quantify in each dataset
combined.peaks <- GenomicRanges::reduce(x = c(gr.PG2, gr.PG6, gr.PG24, gr.PG13))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]



###################################################################################################################################################
#                                                      Creation of seurat objects                                                                  #
###################################################################################################################################################


## Creation PG2 seurat
# create a Seurat object containing the RNA data
soupx_output.PG2 <-readRDS("your/path/to/SoupX/soupx_corrected_countsPG2RN73.rds") #load SoupX contamination corrected output
counts.PG2 <- Read10X_h5("your/path/to/PG2/outs/filtered_feature_bc_matrix.h5") #count data


PG2.dat <- CreateSeuratObject(
  project = "PG2",
  counts = counts.PG2$`Gene Expression`,
  assay = "RNA"
)

## create ATAC assay and add it to the object
#Create Fragment objects
# load metadata
metadata_cellranger.PG2 <- read.table(
  file = "your/path/to/PG2/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

# perform an initial filtering of low count cells
metadata_cellranger.PG2 <- metadata_cellranger.PG2[metadata_cellranger.PG2$is_cell > 0, ]

# create fragment objects
frags.PG2 <- CreateFragmentObject(
  path = "your/path/to/PG2/outs/atac_fragments.tsv.gz",
  cells = rownames(metadata_cellranger.PG2)
)

##quantify peaks
PG2.counts <- FeatureMatrix(
  fragments = frags.PG2,
  features = combined.peaks,
  cells = rownames(metadata_cellranger.PG2)
)

##create object
PG2.dat[["ATAC"]] <- CreateChromatinAssay(
  PG2.counts, 
  fragments = frags.PG2, 
  project = "PG2", 
  assay = "ATAC", 
  sep = c(":", "-"),
  meta.data=metadata_cellranger.PG2,
  annotation = annotation)

#Create corrected RNA data and add to object
PG2.dat[["SoupXRNA"]]<-CreateAssayObject(
  counts=soupx_output.PG2)

#QC cells
DefaultAssay(PG2.dat) <- "ATAC"
PG2.dat <- NucleosomeSignal(PG2.dat)
PG2.dat <- TSSEnrichment(PG2.dat)
PG2.dat<-AddMetaData(PG2.dat,metadata=metadata_cellranger.PG2)

plt.PG2<-VlnPlot(
  object = PG2.dat,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)


saveRDS(PG2.dat,"your/path/to/PG2/outs/PG2.SeuratObject.rds")


## Creation PG24 seurat
# create a Seurat object containing the RNA data
soupx_output.PG24 <-readRDS("your/path/to/SoupX/soupx_corrected_countsPG24RN73.rds") #load SoupX contamination corrected output
counts.PG24 <- Read10X_h5("your/path/to/PG22/outs/filtered_feature_bc_matrix.h5") #count data
PG24.dat <- CreateSeuratObject(
  project = "PG24",
  counts = counts.PG24$`Gene Expression`,
  assay = "RNA"
)

## create ATAC assay and add it to the object
#Create Fragment objects
# load metadata
metadata_cellranger.PG24 <- read.table(
  file = "your/path/to/PG24/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

# perform an initial filtering of low count cells
metadata_cellranger.PG24 <- metadata_cellranger.PG24[metadata_cellranger.PG24$is_cell > 0, ]

# create fragment objects
frags.PG24 <- CreateFragmentObject(
  path = "your/path/to/PG24/outs/atac_fragments.tsv.gz",
  cells = rownames(metadata_cellranger.PG24)
)

##quantify peaks
PG24.counts <- FeatureMatrix(
  fragments = frags.PG24,
  features = combined.peaks,
  cells = rownames(metadata_cellranger.PG24)
)

##create object
PG24.dat[["ATAC"]] <- CreateChromatinAssay(
  PG24.counts, 
  fragments = frags.PG24, 
  project = "PG24", 
  assay = "ATAC", 
  sep = c(":", "-"),
  meta.data=metadata_cellranger.PG24,
  annotation = annotation)

#Create corrected RNA data and add to object
PG24.dat[["SoupXRNA"]]<-CreateAssayObject(
  counts=soupx_output.PG24)

#QC cells
DefaultAssay(PG24.dat) <- "ATAC"
PG24.dat <- NucleosomeSignal(PG24.dat)
PG24.dat <- TSSEnrichment(PG24.dat)
PG24.dat<-AddMetaData(PG24.dat,metadata=metadata_cellranger.PG24)

plt.PG24<-VlnPlot(
  object = PG24.dat,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

saveRDS(PG24.dat,"your/path/to/PG24/outs/PG24.SeuratObject.rds")


## Creation PG6 seurat
# create a Seurat object containing the RNA data
soupx_output.PG6 <-readRDS("your/path/to/SoupX/soupx_corrected_countsPG6RN73.rds") #load SoupX contamination corrected output
counts.PG6 <- Read10X_h5("your/path/to/PG6/outs/filtered_feature_bc_matrix.h5") #count data
PG6.dat <- CreateSeuratObject(
  project = "PG6",
  counts = counts.PG6$`Gene Expression`,
  assay = "RNA"
)

## create ATAC assay and add it to the object
#Create Fragment objects
# load metadata
metadata_cellranger.PG6 <- read.table(
  file = "your/path/to/PG6/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

# perform an initial filtering of low count cells
metadata_cellranger.PG6 <- metadata_cellranger.PG6[metadata_cellranger.PG6$is_cell > 0, ]

# create fragment objects
frags.PG6 <- CreateFragmentObject(
  path = "your/path/to/PG6/outs/atac_fragments.tsv.gz",
  cells = rownames(metadata_cellranger.PG6)
)

##quantify peaks
PG6.counts <- FeatureMatrix(
  fragments = frags.PG6,
  features = combined.peaks,
  cells = rownames(metadata_cellranger.PG6)
)

##create object
PG6.dat[["ATAC"]] <- CreateChromatinAssay(
  PG6.counts, 
  fragments = frags.PG6, 
  project = "PG6", 
  assay = "ATAC", 
  sep = c(":", "-"),
  meta.data=metadata_cellranger.PG6,
  annotation = annotation)

#Create corrected RNA data and add to object
PG6.dat[["SoupXRNA"]]<-CreateAssayObject(
  counts=soupx_output.PG6)

#QC cells
DefaultAssay(PG6.dat) <- "ATAC"
PG6.dat <- NucleosomeSignal(PG6.dat)
PG6.dat <- TSSEnrichment(PG6.dat)
PG6.dat<-AddMetaData(PG6.dat,metadata=metadata_cellranger.PG6)

plt.PG6<-VlnPlot(
  object = PG6.dat,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

saveRDS(PG6.dat,"your/path/to/PG6/outs/PG6.SeuratObject.rds")

## Creation PG13 seurat
# create a Seurat object containing the RNA data
soupx_output.PG13 <-readRDS("your/path/to/SoupX/soupx_corrected_countsPG13RN73.rds") 
counts.PG13 <- Read10X_h5("your/path/to/PG13/outs/filtered_feature_bc_matrix.h5")
PG13.dat <- CreateSeuratObject(
  project = "PG13",
  counts = counts.PG13$`Gene Expression`,
  assay = "RNA"
)

## create ATAC assay and add it to the object
#Create Fragment objects
# load metadata
metadata_cellranger.PG13 <- read.table(
  file = "your/path/to/PG13/outs/per_barcode_metrics.csv",
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1, ] # remove the first row

# perform an initial filtering of low count cells
metadata_cellranger.PG13 <- metadata_cellranger.PG13[metadata_cellranger.PG13$is_cell > 0, ]

# create fragment objects
frags.PG13 <- CreateFragmentObject(
  path = "your/path/to/PG13/outs/atac_fragments.tsv.gz",
  cells = rownames(metadata_cellranger.PG13)
)

##quantify peaks
PG13.counts <- FeatureMatrix(
  fragments = frags.PG13,
  features = combined.peaks,
  cells = rownames(metadata_cellranger.PG13)
)

##create object
PG13.dat[["ATAC"]] <- CreateChromatinAssay(
  PG13.counts, 
  fragments = frags.PG13, 
  project = "PG13", 
  assay = "ATAC", 
  sep = c(":", "-"),
  meta.data=metadata_cellranger.PG13,
  annotation = annotation)

#Create corrected RNA data and add to object
PG13.dat[["SoupXRNA"]]<-CreateAssayObject(
  counts=soupx_output.PG13)

#QC cells
DefaultAssay(PG13.dat) <- "ATAC"
PG13.dat <- NucleosomeSignal(PG13.dat)
PG13.dat <- TSSEnrichment(PG13.dat)
PG13.dat <- AddMetaData(PG13.dat,metadata=metadata_cellranger.PG13)


plt.PG13<-VlnPlot(
  object = PG13.dat,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)

saveRDS(PG13.dat,"your/path/to/PG13/outs/PG13.SeuratObject.rds")


#  %MT Filtre & RNA count

PG2.dat <- readRDS("your/path/to/PG2/outs/PG2.SeuratObject.rds")
PG6.dat <- readRDS("your/path/to/PG6/outs/PG6.SeuratObject.rds")
PG24.dat <- readRDS("your/path/to/PG24/outs/PG24.SeuratObject.rds")
PG13.dat <- readRDS("your/path/to/PG13/outs/PG13.SeuratObject.rds")


DefaultAssay(PG2.dat) <- "SoupXRNA"
DefaultAssay(PG6.dat) <- "SoupXRNA"
DefaultAssay(PG24.dat) <- "SoupXRNA"
DefaultAssay(PG13.dat) <- "SoupXRNA"


mt_genes <- c("ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3", 
              "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB")


PG2.dat[["percent.mtPG2"]] <-PercentageFeatureSet(PG2.dat, features = mt_genes)
PG2.dat <- subset(PG2.dat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mtPG2 < 1)

PG6.dat[["percent.mtPG6"]] <-PercentageFeatureSet(PG6.dat, features = mt_genes)
PG6.dat <- subset(PG6.dat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mtPG6 < 1)

PG24.dat[["percent.mtPG24"]] <-PercentageFeatureSet(PG24.dat, features = mt_genes)
PG24.dat <- subset(PG24.dat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mtPG24 < 1)

PG13.dat[["percent.mtPG13"]] <-PercentageFeatureSet(PG13.dat, features = mt_genes)
PG13.dat <- subset(PG13.dat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mtPG13 < 1)


# Computation of LSI for ATAC-seq Data 
DefaultAssay(PG2.dat) <- "ATAC"
DefaultAssay(PG6.dat) <- "ATAC"
DefaultAssay(PG24.dat) <- "ATAC"
DefaultAssay(PG13.dat) <- "ATAC"

PG2.dat <- FindTopFeatures(PG2.dat, min.cutoff = 10)
PG2.dat <- RunTFIDF(PG2.dat)
PG2.dat <- RunSVD(PG2.dat)

PG6.dat <- FindTopFeatures(PG6.dat, min.cutoff = 10)
PG6.dat <- RunTFIDF(PG6.dat)
PG6.dat <- RunSVD(PG6.dat)

PG24.dat <- FindTopFeatures(PG24.dat, min.cutoff = 10)
PG24.dat <- RunTFIDF(PG24.dat)
PG24.dat <- RunSVD(PG24.dat)

PG13.dat <- FindTopFeatures(PG13.dat, min.cutoff = 10)
PG13.dat <- RunTFIDF(PG13.dat)
PG13.dat <- RunSVD(PG13.dat)


###################################################################################################################################################
#                                                      Sample integration                                                                         #
###################################################################################################################################################


## Merge all samples in an unic object
# rename sample
PG2.dat <- RenameCells(
  PG2.dat,
  new.names = paste0("D2_", colnames(x =PG2.dat))
)

PG6.dat <- RenameCells(
  PG6.dat,
  new.names = paste0("PM_", colnames(x =PG6.dat))
)

PG24.dat <- RenameCells(
  PG24.dat,
  new.names = paste0("PS_", colnames(x =PG24.dat))
)

PG13.dat <- RenameCells(
  PG13.dat,
  new.names = paste0("E_", colnames(x =PG13.dat))
)

# add information to identify dataset of origin
PG2.dat$dataset <- 'D2'
PG6.dat$dataset <- 'PM'
PG24.dat$dataset <- 'PS'
PG13.dat$dataset <- 'E'

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = PG2.dat,
  y = list(PG6.dat, PG24.dat, PG13.dat),
)
combined[["ATAC"]]

combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')

P1 <- DimPlot(combined, group.by = 'dataset', pt.size = 0.1)


### ATAC integration 
# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = list(PG2.dat, PG6.dat, PG24.dat, PG13.dat),
  anchor.features = rownames(D2),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 2:30
)

# create a new UMAP using the integrated embeddings
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:30)
p2 <- DimPlot(integrated73, group.by = "dataset")
p2 + ggtitle("integrated73")
(P1 + ggtitle("Merged")) | (p2 + ggtitle("integrated"))

saveRDS(integrated, file='your/path/to/PGintegrated.rds')




### Call Peaks and Dimensionality Reduction

# call peaks using MACS3
DefaultAssay(integrated)<-"ATAC"
peaks.int <- CallPeaks(integrated, macs2.path = 'your/path/to/macs3')


# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(integrated),
  features = peaks.int,
  cells = colnames(integrated)
)

# create a new assay using the MACS3 peak set and add it to the Seurat object
### get gene annotations for Rn7
integrated[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = integrated@assays$ATAC@fragments,
  annotation = annotation
)


#set up colors for samples
my_cols = brewer.pal(1,"Spectral")
alpha_val=0.33


#RNA Processing
DefaultAssay(integrated) <- "SoupXRNA"

# Run the standard workflow for visualization and clustering
# Filtering
integrated[["percent.mt"]] <-PercentageFeatureSet(integrated,features = mt_genes)                        
integrated <- subset(integrated73, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 0.1)

#Normalisarion and dimension reduction
integrated <- NormalizeData(integrated)
integrated <- SCTransform(integrated, vars.to.regress = "percent.mt", verbose = FALSE)
integrated <- FindVariableFeatures(integrated)
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated)
ElbowPlot(integrated)

# UMAP
integrated <- RunUMAP(
  object = integrated,
  reduction.name="RNA_umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=2:30
)

p1a <-DimPlot(integrated,reduction="RNA_umap")+ggtitle("RNA RNA_umap UMAP")
p1b <-DimPlot(integrated,reduction="RNA_umap", group.by="dataset")

integrated$sample_Id <- integrated$orig.ident
integrated$orig.ident <- integrated$dataset

# First gene visualisation

DefaultAssay(integrated) <- "SoupXRNA"

#MT
FeaturePlot(integrated, reduction="RNA_umap", features = c( "COX2", "ATP6", "ND1"), min.cutoff = 'q10')
#Somatotropes
FeaturePlot(integrated73, reduction="RNA_umap", features = c( "Gh1", "Prl", "S100g", "Pou1f1", "Tpt1", "Gal", "Gria4"), max.cutoff = 'q90', order = TRUE)
#Lactotropes
FeaturePlot(integrated73, reduction="RNA_umap", features = c( "Prl", "S100g"), order = TRUE)
#Gonadotropes
FeaturePlot(integrated73, reduction="RNA_umap", features = c( "Lhb", "Fshb", "Nr5a1"), order = TRUE)


#DNA Accessibility processing
DefaultAssay(integrated) <- "peaks"
integrated <- FindTopFeatures(integrated, min.cutoff = 5)
integrated <- RunTFIDF(integrated)
integrated <- RunSVD(integrated)
integrated<- RunUMAP(
  object = integrated,
  reduction.name="atac_umap",
  reduction="lsi",
  assay = "peaks",
  verbose = TRUE,
  dims=2:40
)
p2<-DimPlot(integrated,reduction="atac_umap")+ggtitle("ATAC UMAP")



# build a multimodal graph
integrated <- FindMultiModalNeighbors(
  object = integrated,
  reduction.list = list("pca", "lsi" ), 
  dims.list = list(3:15, 1:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
integrated <- RunUMAP(
  object = integrated,
  nn.name = "weighted.nn",
  reduction.name="multimodal_umap",
  assay = "RNA",
  verbose = TRUE,
  return.neighbor = TRUE
)

p3<-DimPlot(integrated,reduction="multimodal_umap",group.by="dataset")


#Cluster on multimodal graph
integrated <- FindClusters(integrated, resolution = 0.4, verbose = FALSE, graph="wknn")
p4<-DimPlot(integrated,reduction="multimodal_umap", group.by = "seurat_clusters", label = TRUE)
p5<-DimPlot(integrated,reduction="RNA_umap", group.by = "seurat_clusters", label = TRUE)
p6<-DimPlot(integrated,reduction="atac_umap", group.by = "seurat_clusters", label = TRUE)


#Cluster identification
DefaultAssay(integrated) <- "SoupXRNA"

# Markers of Somatotrope and Lactotrope cells
FeaturePlot(integrated, reduction="multimodal_umap", features = c("Gh1", "Prl", "Ghrhr", "Pappa2", "Agtr1b", "Alk"), max.cutoff = 'q90')
# Markers of Gonadotrope cells
FeaturePlot(integrated, reduction="multimodal_umap", features = c("Lhb", "Fshb", "Cga", "Gnrhr","Nr5a1"), max.cutoff = 'q90')
# Markers of Thyrotrope cells
FeaturePlot(integrated, reduction="multimodal_umap", features = c("Tshb", "Cga"), max.cutoff = 'q90')
# Markers of Corticotrope and Melanotrope cells
FeaturePlot(integrated, reduction="multimodal_umap", features = c("Pomc", "Pax7", "Neurod1", "Ly6h"), max.cutoff = 'q90')
# Markers of Folliculostellate cells
FeaturePlot(integrated, reduction="multimodal_umap", features = c("Adamts9", "Rfx4", "S100b", "Prop1"), max.cutoff = 'q90')
# Markers of Endothelial cells
FeaturePlot(integrated, reduction="multimodal_umap", features = c("Plvap", "Emcn"), max.cutoff = 'q90')
# Markers of Leucocytes
FeaturePlot(integrated, reduction="multimodal_umap", features = c("Tyrobp", "Arhgap15"), max.cutoff = 'q90')
# Markers of Pericytes
FeaturePlot(integrated, reduction="multimodal_umap", features = c("Pde5a", "Pdgfrb"), max.cutoff = 'q90')
# Markers of Pericytes
FeaturePlot(integrated, reduction="multimodal_umap", features = c("Lmod3"), max.cutoff = 'q90')


#### remove low quality cluster
integrated <- subset(integrated, idents = 6, invert = TRUE)

# CLUSTERS identification
# name the clusters
integrated <- RenameIdents(integrated, 
                           "0" = "Cluster_S",
                           "1" = "Cluster_L", 
                           "2" = "Cluster_L",
                           "3" = "Cluster_FSC",  
                           "4" = "Cluster_L",
                           "5" = "Cluster_C",
                           "7" = "Cluster_Le", 
                           "8" = "Cluster_G", 
                           "9" = "Cluster_T",
                           "10" = "Cluster_M",
                           "11" = "Cluster_G", 
                           "12" = "Cluster_EC",
                           "13" = "Cluster_S", 
                           "14" = "Cluster_Pe"
)

# add cell_type metadata
cluster_names <- Idents(integrated)
integrated <- AddMetaData(
  object =  integrated,
  metadata = cluster_names,
  col.name = 'cell_Type')


# Extract SoupX corrected RNA counts
DefaultAssay(integrated) <- "SoupXRNA"
countsSoupXRNA <- AggregateExpression(integrated, group.by = c("cell_Type", "dataset"),
                                      assays = 'SoupXRNA',
                                      slot = "counts"
)
write.csv(countsSoupXRNA,"your/path/to/counts_integrated.csv")


# Save final objet
save.rds(integrated, "your/path/to/PGintegrated.rds")

# Create gene expression and activity plot
GeneOfInterest <- "Prl"
DefaultAssay(integrated) <- "SoupXRNA"
FeaturePlot_RNA_GeneOfInterest <- FeaturePlot(integrated, reduction="multimodal_umap", features = GeneOfInterest, cols = c("#eeeeee", "#90255a"), pt.size = 1
                                              , order = TRUE , min.cutoff = 'q10'   , max.cutoff = 'q90')
ggsave('./FeaturePlot_RNA_GeneOfInterest.png', FeaturePlot_RNA_GeneOfInterest)
DefaultAssay(integrated) <- "GeneActivity"
FeaturePlot_GA_GeneOfInterest <- FeaturePlot(integrated, reduction="multimodal_umap", features = GeneOfInterest, cols = c("#eeeeee", "#90255a"), pt.size = 1
                                             , order = TRUE , min.cutoff = 'q10'   , max.cutoff = 'q90')
ggsave('./FeaturePlot_GA_GeneOfInterest.png', FeaturePlot_GA_GeneOfInterest)


###################################################################################################################################################
#                                                     Link peaks to genes                                                                         #
###################################################################################################################################################


library(BSgenome)
library(BSgenome.Rnorvegicus.NCBI.rn7.2)
genome <- BSgenome.Rnorvegicus.NCBI.rn7.2
seqlevelsStyle(genome) <- "UCSC"


integrated$PG_cell_Type <- paste0(integrated$dataset, "_", integrated$cell_Type)
Idents(integrated) <- "dataset"



##### Load fragment tables
integrated@assays$peaks@fragments[[1]]@path <- 'your/path/to/PG2/outs/atac_fragments.tsv.gz'
integrated@assays$peaks@fragments[[2]]@path <- 'your/path/to/PG6/outs/atac_fragments.tsv.gz'
integrated@assays$peaks@fragments[[3]]@path <- 'your/path/to/PG24/outs/atac_fragments.tsv.gz'
integrated@assays$peaks@fragments[[4]]@path <- 'your/path/to/PG13/outs/atac_fragments.tsv.gz'

##Linking peaks to genes
DefaultAssay(integrated)<-"peaks"
main.chroms <- standardChromosomes(BSgenome.Rnorvegicus.NCBI.rn7.2)
keep.peaks <- as.logical(seqnames(granges(integrated)) %in% main.chroms)
integrated[["peaks"]] <- subset(integrated[["peaks"]], features = rownames(integrated[["peaks"]])[keep.peaks])
integrated <- RegionStats(integrated, genome = BSgenome.Rnorvegicus.NCBI.rn7.2)

# link peaks to genes
DefaultAssay(integrated)<-"peaks"
integrated <- LinkPeaks(
  integrated,
  peak.assay = 'peaks',
  expression.assay  = "SoupXRNA",
  peak.slot = "counts",
  expression.slot = "data",
  method = "pearson",
  gene.coords = NULL,
  distance = 5.2e+06,
  min.distance = NULL,
  min.cells = 10,
  genes.use = NULL,
  n_sample = 1000,
  pvalue_cutoff = 0.05,
  score_cutoff = 0.05,
  gene.id = FALSE,
  verbose = TRUE
)

write.csv(integrated@assays$peaks$links,"your/path/to/linkpeaks_integrated_all_cluster.csv")


save.rds(integrated, "your/path/to/PGintegrated.rds")