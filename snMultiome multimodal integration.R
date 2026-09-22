
################################################################################
## 00. INITIALISATION
################################################################################

rm(list = ls())
gc()

set.seed(1234)

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({

  library(Seurat)
  library(Signac)

  library(GenomicRanges)
  library(GenomeInfoDb)
  library(EnsDb.Rnorvegicus.v110)

  library(SoupX)

  library(ggplot2)
  library(patchwork)

  library(data.table)
  library(Matrix)

  library(future)
  library(harmony)
})

future::plan(
    "multisession",
    workers = 8
)

options(
    future.globals.maxSize = 30000 * 1024^2
)

options(future.seed = TRUE)

################################################################################
## 01. PATHS
################################################################################

project.dir <- "/shared/projects/femcycle/PG_multiome_cluster"

data.dir <- file.path(project.dir, "data")

results.dir <- file.path(project.dir, "results")

figures.dir <- file.path(results.dir, "comparison73", "figures")

tables.dir <- file.path(results.dir, "comparison73", "tables")

objects.dir <- file.path(results.dir, "comparison73", "objects")

qc.dir <- file.path(results.dir, "comparison73", "qc")

dir.create(figures.dir, recursive = TRUE, showWarnings = FALSE)

dir.create(tables.dir, recursive = TRUE, showWarnings = FALSE)

dir.create(objects.dir, recursive = TRUE, showWarnings = FALSE)

dir.create(qc.dir, recursive = TRUE, showWarnings = FALSE)

dir.create(
    file.path(qc.dir, "soupx"),
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    file.path(qc.dir, "cell_filter"),
    recursive = TRUE,
    showWarnings = FALSE
)

macs2.path <- Sys.which("macs3")

if (macs2.path == "")
    stop("macs3 not found in PATH")

################################################################################
## 02. PARAMETERS
################################################################################

samples <- c(
  "PG1",
  "PG2",
  "PG5",
  "PG6",
  "PG9",
  "PG11",
  "PG13",
  "PG24",
  "PG32",
  "PG38",
  "PGo1",
  "PGo6"
)

sample.condition <- c(
    PG1  = "E",
    PG2  = "D2",
    PG5  = "PM",
    PG6  = "PM",
    PG9  = "E",
    PG11 = "E",
    PG13 = "E",
    PG24 = "PS",
    PG32 = "PS",
    PG38 = "PS",
    PGo1 = "D2",
    PGo6 = "D2"
)

cellranger.dir <- c(
  PG1  = file.path(project.dir, "PG1RN73"),
  PG2  = file.path(project.dir, "PG2RN73"),
  PG5  = file.path(project.dir, "PG5RN73"),
  PG6  = file.path(project.dir, "PG6RN73"),
  PG9  = file.path(project.dir, "PG9RN73"),
  PG11 = file.path(project.dir, "PG11RN73"),
  PG13 = file.path(project.dir, "PG13RN73"),
  PG24 = file.path(project.dir, "PG24RN73"),
  PG32 = file.path(project.dir, "PG32RN73"),
  PG38 = file.path(project.dir, "PG38RN73"),
  PGo1 = file.path(project.dir, "PGo1RN73"),
  PGo6 = file.path(project.dir, "PGo6RN73")
)


for (s in samples) {

    outs <- file.path(cellranger.dir[[s]], "outs")

    required <- c(
        "atac_fragments.tsv.gz",
        "atac_fragments.tsv.gz.tbi",
        "raw_feature_bc_matrix.h5",
        "filtered_feature_bc_matrix.h5",
        "per_barcode_metrics.csv"
    )

    missing <- required[
        !file.exists(file.path(outs, required))
    ]

    if (length(missing) > 0) {
        stop(
            s,
            " : missing files: ",
            paste(missing, collapse = ", ")
        )
    }
}



call_macs3_sample <- function(sample, sample_dir) {

    outs <- file.path(
        sample_dir,
        "outs"
    )

    frag <- file.path(
        outs,
        "atac_fragments.tsv.gz"
    )

    tmp.dir <- file.path(
        results.dir,
        "peak_calling_by_sample",
        "fragment_tmp"
    )

    macs.dir <- file.path(
        results.dir,
        "peak_calling_by_sample",
        "macs_output"
    )

    dir.create(
        tmp.dir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    dir.create(
        macs.dir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    bed <- file.path(
        tmp.dir,
        paste0(sample, ".bed")
    )

    message("\n================================")
    message("MACS3 : ", sample)
    message("================================")

    ############################################################
    ## Fragment file -> BED
    ############################################################

    cmd <- paste(
        "zcat",
        shQuote(frag),
        "| awk 'BEGIN{OFS=\"\\t\"} {print $1,$2,$3}'",
        ">",
        shQuote(bed)
    )

    status <- system(cmd)

    if (
        status != 0 ||
        !file.exists(bed) ||
        file.size(bed) == 0
    ) {
        stop(
            "Failed to create BED for ",
            sample
        )
    }

    ############################################################
    ## MACS3
    ############################################################

    args <- c(
        "callpeak",
        "-t", bed,
        "-g", "2.7e9",
        "-f", "BED",
        "--nomodel",
        "--extsize", "200",
        "--shift", "-100",
        "-n", sample,
        "--outdir", macs.dir
    )

    status <- system2(
        macs2.path,
        args = args
    )

    if (status != 0) {
        stop(
            "MACS3 failed for ",
            sample
        )
    }

    peak.file <- file.path(
        macs.dir,
        paste0(sample, "_peaks.narrowPeak")
    )

    if (
        !file.exists(peak.file) ||
        file.size(peak.file) == 0
    ) {
        stop(
            "MACS3 output missing for ",
            sample,
            ": ",
            peak.file
        )
    }


    peaks.dt <- data.table::fread(
        peak.file,
        header = FALSE
    )

    if (ncol(peaks.dt) < 3) {
        stop(
            "Invalid MACS3 narrowPeak file for ",
            sample
        )
    }

    peaks <- GenomicRanges::GRanges(
        seqnames = peaks.dt[[1]],
        ranges = IRanges::IRanges(
            start = as.integer(peaks.dt[[2]]) + 1L,
            end = as.integer(peaks.dt[[3]])
        ),
        strand = "*"
    )

    if (ncol(peaks.dt) >= 10) {

        mcols(peaks)$name <-
            as.character(peaks.dt[[4]])

        mcols(peaks)$score <-
            as.numeric(peaks.dt[[5]])

        mcols(peaks)$signalValue <-
            as.numeric(peaks.dt[[7]])

        mcols(peaks)$pValue <-
            as.numeric(peaks.dt[[8]])

        mcols(peaks)$qValue <-
            as.numeric(peaks.dt[[9]])

        mcols(peaks)$peak <-
            as.integer(peaks.dt[[10]])

    }

    rm(peaks.dt)

    message(
        sample,
        " : ",
        length(peaks),
        " peaks"
    )

    peaks
}


soupx.contamination <- c(

  PG2  = 0.15,

  PG6  = 0.15,

  PG13 = 0.08,

  PG24 = 0.15,

  PG1 = 0.1,
  PG5 = 0.25,
  PG9 = 0.15,
  PG11 = 0.3,
  PG32 = 0.2,
  PG38 = 0.18,
  PGo1 = 0.15,
  PGo6 = 0.20

)

mt.threshold <- c(

    PG2  = 1,

    PG6  = 1,

    PG13 = 1,

    PG24 = 1,

    PG1 = 3,
    PG5 = 3,
    PG9 = 3,
    PG11 = 3,
    PG32 = 3,
    PG38 = 3,
    PGo1 = 3,
    PGo6 = 3

)

rna.min.features <- 200

rna.max.features <- 4500

lsi.min.cutoff <- 10

soupx.markers <- c(

    "Lhb",
    "Fshb",
    "Prl",
    "Pomc",
"Nr5a1",
"Gnrhr",
"S100g",
"Cga",
"S100b",
"Tshb",
"Tyrobp",
"Plvap"

)

reference <- file.path(data.dir)

old.object <- file.path(

  results.dir,

  "PGintegrated73.reclustered.rds"

)

new.object <- file.path(

  objects.dir,

  "PGintegrated73_v2.rds"

)

cluster.annotation <- c(
#Cluster_S   Cluster_L Cluster_FSC   Cluster_C  Cluster_Le   Cluster_G   Cluster_T   Cluster_M  Cluster_EC  Cluster_Pe
    "0"  = "Cluster_S",
    "1"  = "Cluster_L",
    "2"  = "Cluster_L",
    "3"  = "Cluster_FSC",
    "4"  = "Cluster_L",
    "5"  = "Cluster_S",
    "6"  = "Cluster_L",
    "7"  = "Cluster_C",
    "8"  = "Cluster_G",
    "9"  = "Cluster_Le",
    "10" = "Cluster_L",
    "11" = "Cluster_M",
    "12" = "Cluster_G",
    "13" = "Cluster_L",
    "14" = "Cluster_EC",
    "15" = "Cluster_T",
    "16" = "Cluster_Pe"
)


################################################################################
## 03. FUNCTIONS
################################################################################

read_cellbender_h5 <- function(file){

    f <- hdf5r::H5File$new(
        file,
        mode = "r"
    )

    g <- f[["matrix"]]

    mat <- Matrix::sparseMatrix(

        i = g[["indices"]]$read() + 1,

        p = g[["indptr"]]$read(),

        x = g[["data"]]$read(),

        dims = g[["shape"]]$read()

    )

    rownames(mat) <- g[["features"]][["name"]]$read()

    colnames(mat) <- g[["barcodes"]]$read()

    f$close_all()

    mat
}


load_sample <- function(sample){

  sample.path <- cellranger.dir[[sample]]

  if (is.null(sample.path))
    stop("Unknown sample: ", sample)

  outs <- file.path(sample.path, "outs")

  message("Loading ", sample)

  ## RNA

  ## Cell Ranger matrices

  raw <- Read10X_h5(
    file.path(
      outs,
      "raw_feature_bc_matrix.h5"
    )
  )

  filtered <- Read10X_h5(
    file.path(
      outs,
      "filtered_feature_bc_matrix.h5"
    )
  )

  rna.raw <- raw[["Gene Expression"]]

  atac.raw <- raw[["Peaks"]]

atac.filtered <- filtered[["Peaks"]]



###################

rna.filtered <- filtered[["Gene Expression"]]


##################


################################################################################
## CELLBENDER FPR 0.001 + CELL RANGER ARC CELL FILTER
################################################################################

cellbender.file <- file.path(
    outs,
    paste0(
        "cellbender_",
        sample,
        "_FPR_0.001_filtered.h5"
    )
)

if(!file.exists(cellbender.file)){

    stop(
        "CellBender FPR 0.001 file not found for ",
        sample,
        ": ",
        cellbender.file
    )

}

cellbender.rna <- read_cellbender_h5(
    cellbender.file
)

cellranger.cells <- colnames(
    filtered[["Gene Expression"]]
)

common.cells <- intersect(
    cellranger.cells,
    colnames(cellbender.rna)
)

message(
    sample,
    " : Cell Ranger ARC cells = ",
    length(cellranger.cells),
    " ; CellBender cells = ",
    ncol(cellbender.rna),
    " ; common = ",
    length(common.cells)
)

if(length(common.cells) == 0){

    stop(
        "No common cells between Cell Ranger ARC and CellBender for ",
        sample
    )

}

rna.filtered <- cellbender.rna[
    ,
    common.cells,
    drop = FALSE
]

################################################################################
## ALIGN GENES RAW / CELLBENDER
################################################################################

common.genes <- intersect(
    rownames(rna.raw),
    rownames(rna.filtered)
)

rna.raw <- rna.raw[
    common.genes,
    ,
    drop = FALSE
]

rna.filtered <- rna.filtered[
    common.genes,
    ,
    drop = FALSE
]

rna.filtered <- rna.filtered[
    rownames(rna.raw),
    ,
    drop = FALSE
]

#####################################





  ## ATAC peaks

  peaks <- rtracklayer::import(
    file.path(
      outs,
      "atac_peaks.bed"
    )
  )

  ## Metrics

  metrics <- data.table::fread(
    file.path(
        outs,
        "per_barcode_metrics.csv"
    )
)

metrics <- as.data.frame(metrics)

rownames(metrics) <- metrics[[1]]

metrics <- metrics[, -1, drop = FALSE]

  summary <- data.table::fread(
    file.path(
      outs,
      "summary.csv"
    )
  )


stopifnot(
    identical(
        rownames(rna.raw),
        rownames(rna.filtered)
    )
)

message(
    sample,
    " : final RNA = ",
    nrow(rna.filtered),
    " genes x ",
    ncol(rna.filtered),
    " cells"
)




   return(list(

    name = sample,

    path = sample.path,

    rna.raw = rna.raw,

    rna.filtered = rna.filtered,

    atac.raw = atac.raw,

    atac.filtered = atac.filtered,

    fragments = file.path(
      outs,
      "atac_fragments.tsv.gz"
    ),

    peaks = peaks,

    metrics = metrics,

    summary = summary

  ))

}



run_soupx <- function(sample){

    message("Running SoupX: ", sample$name)

    srat <- CreateSeuratObject(
        counts = sample$rna.filtered
    )

    soup <- SoupChannel(
        sample$rna.raw,
        sample$rna.filtered
    )

oplan <- future::plan()
future::plan("sequential")

srat <- SCTransform(
    srat,
    verbose = FALSE
)

future::plan(oplan)

    srat <- RunPCA(
        srat,
        verbose = FALSE
    )

    srat <- RunUMAP(
        srat,
        dims = 1:15,
        verbose = FALSE
    )

    srat <- FindNeighbors(
        srat,
        dims = 1:15,
        verbose = FALSE
    )

    srat <- FindClusters(
        srat,
        verbose = FALSE
    )

    soup <- setClusters(
        soup,
        srat$seurat_clusters
    )

    soup <- setDR(
        soup,
        Embeddings(srat, "umap")
    )

cat("\n==============================\n")
cat("SoupX diagnostics :", sample$name, "\n")
cat("==============================\n")
cat("Raw cells      :", ncol(sample$rna.raw), "\n")
cat("Filtered cells :", ncol(sample$rna.filtered), "\n")
cat("Genes          :", nrow(sample$rna.filtered), "\n")
cat("Clusters       :", length(unique(srat$seurat_clusters)), "\n")
cat("Cells/cluster\n")
print(table(srat$seurat_clusters))
cat("==============================\n")

print(head(rownames(sample$rna.filtered), 20))

   message("Running SoupX autoEstCont (QC only)")

tryCatch({

    soup <- autoEstCont(soup)

}, error = function(e) {

    warning(
        sample$name,
        ": autoEstCont failed: ",
        conditionMessage(e)
    )

    soup$fit <- NULL

})

message("Continuing with manual contamination fraction")

    sample$soupx <- soup

    sample$quick.umap <- srat

    return(sample)
}



plot_rho_qc <- function(sample){

    message("SoupX rho: ", sample$name)

    fit <- sample$soupx$fit

if (is.null(fit) || is.null(fit$rhoEst)) {

    cat(
        "\n=============================\n",
        "Sample: ", sample$name, "\n",
        "autoEstCont failed\n",
        "Applied rho   : ", soupx.contamination[[sample$name]], "\n",
        "=============================\n"
    )

    return(invisible(sample))
}

    cat(
        "\n=============================\n",
        "Sample: ", sample$name, "\n",
        "Estimated rho : ", signif(fit$rhoEst, 4), "\n",
        "Applied rho   : ", soupx.contamination[[sample$name]], "\n",
        "Prior rho     : ", fit$priorRho, "\n",
        "95% interval  : ",
        signif(fit$rhoFWHM[1], 4),
        " - ",
        signif(fit$rhoFWHM[2], 4),
        "\n",
        "Markers used  : ", nrow(fit$markersUsed), "\n",
        "=============================\n"
    )

    write.table(

        data.frame(

            sample = sample$name,

            estimated_rho = fit$rhoEst,

            applied_rho = soupx.contamination[[sample$name]],

            rho_low = fit$rhoFWHM[1],

            rho_high = fit$rhoFWHM[2],

            markers_used = nrow(fit$markersUsed)

        ),

        file = file.path(
            qc.dir,
            "soupx",
            paste0(sample$name, "_rho.tsv")
        ),

        sep = "\t",

        row.names = FALSE,

        quote = FALSE

    )

    invisible(sample)

}



apply_soupx <- function(sample){

    rho <- soupx.contamination[[sample$name]]

    if(is.na(rho))
        stop("SoupX contamination not defined for ", sample$name)

    sample$soupx <- setContaminationFraction(

        sample$soupx,

        rho,

        forceAccept = TRUE

    )

    sample$rna.corrected <- adjustCounts(

        sample$soupx,

        roundToInt = TRUE

    )

sample$quick.umap[["SoupXRNA"]] <- CreateAssayObject(
    counts = sample$rna.corrected
)

    sample

}



plot_soupx_qc <- function(sample){

    message("SoupX QC: ", sample$name)

    obj <- sample$quick.umap

    DefaultAssay(obj) <- "SoupXRNA"

    dir.create(
        file.path(qc.dir, "soupx"),
        recursive = TRUE,
        showWarnings = FALSE
    )

    pdf(
        file.path(
            qc.dir,
            "soupx",
            paste0(sample$name, "_SoupX_QC.pdf")
        ),
        width = 10,
        height = 8
    )

markers.per.page <- split(
    soupx.markers,
    ceiling(seq_along(soupx.markers) / 4)
)

for (feat in markers.per.page) {

    print(
        FeaturePlot(
            obj,
            reduction = "umap",
            features = feat,
            order = TRUE,
            ncol = 2,
            max.cutoff = "q95"
        )
    )

}

    dev.off()

    invisible(sample)

}



create_seurat <- function(
    sample,
    combined.peaks,
    annotation
){

    message("Creating Seurat object: ", sample$name)

    metadata <- sample$metrics

    metadata <- metadata[metadata$is_cell > 0, , drop = FALSE]

#stopifnot(
#    all(rownames(metadata) %in% colnames(sample$rna.filtered))
#)

############


################################################################################
## ALIGN CELL RANGER METADATA WITH CELLBENDER CELLS
################################################################################

metadata <- metadata[
    rownames(metadata) %in% colnames(sample$rna.filtered),
    ,
    drop = FALSE
]

sample$rna.filtered <- sample$rna.filtered[
    ,
    rownames(metadata),
    drop = FALSE
]

message(
    sample$name,
    " : metadata cells = ",
    nrow(metadata),
    " ; RNA cells = ",
    ncol(sample$rna.filtered)
)

stopifnot(
    identical(
        rownames(metadata),
        colnames(sample$rna.filtered)
    )
)

###############


    frags <- CreateFragmentObject(
        path = sample$fragments,
        cells = rownames(metadata)
    )

    atac.counts <- FeatureMatrix(
        fragments = frags,
        features = combined.peaks,
        cells = rownames(metadata)
    )

stopifnot(
    ncol(atac.counts) == nrow(metadata)
)

    obj <- CreateSeuratObject(
        counts = sample$rna.filtered,
        assay = "RNA",
        project = sample$name
    )

    obj[["ATAC"]] <- CreateChromatinAssay(
        counts = atac.counts,
        fragments = frags,
        project = sample$name,
        assay = "ATAC",
        sep = c(":", "-"),
        meta.data = metadata,
        annotation = annotation
    )

sample$rna.corrected <-
    sample$rna.corrected[
        ,
        colnames(obj)
    ]

    obj[["SoupXRNA"]] <- CreateAssayObject(
        counts = sample$rna.corrected
    )

    DefaultAssay(obj) <- "ATAC"

    obj <- NucleosomeSignal(obj)

obj <- TSSEnrichment(obj)

    obj <- AddMetaData(
        obj,
        metadata = metadata
    )

stopifnot(

    ncol(obj) == ncol(sample$rna.filtered),

    ncol(obj) == ncol(sample$rna.corrected)

)

    sample$object <- obj

    return(sample)
}


plot_cell_qc <- function(sample){

    message("Cell QC: ", sample$name)

    obj <- sample$object

    mt.genes <- c(

        "ND1",
        "ND2",
        "COX1",
        "COX2",
        "ATP8",
        "ATP6",
        "COX3",
        "ND3",
        "ND4L",
        "ND4",
        "ND5",
        "ND6",
        "CYTB"

    )

    DefaultAssay(obj) <- "SoupXRNA"

    obj[["percent.mt"]] <- PercentageFeatureSet(

        obj,

        features = mt.genes

    )

    p1 <- VlnPlot(

        object = obj,

        features = c(

            "nCount_SoupXRNA",

            "nFeature_SoupXRNA",

            "nCount_ATAC",

            "TSS.enrichment",

            "nucleosome_signal",

            "percent.mt"

        ),

        pt.size = 0,

        ncol = 3

    )

    p2 <- FeatureScatter(

        obj,

        feature1 = "nCount_SoupXRNA",

        feature2 = "percent.mt"

    )

    p3 <- FeatureScatter(

        obj,

        feature1 = "nFeature_SoupXRNA",

        feature2 = "percent.mt"

    )

    pdf(

        file.path(

            qc.dir,

            "cell_filter",

            paste0(sample$name, "_cell_QC.pdf")

        ),

        width = 14,

        height = 12

    )

    print(p1)

    print(p2)

    print(p3)

    dev.off()

    sample$object <- obj

    return(sample)

}

filter_cells <- function(sample){

    message("Filtering cells: ", sample$name)

    obj <- sample$object

if(!"percent.mt" %in% colnames(obj@meta.data))
    stop("Run plot_cell_qc() before filter_cells().")

    DefaultAssay(obj) <- "SoupXRNA"

    obj <- subset(

        obj,

        subset =

            nFeature_SoupXRNA > rna.min.features &

            nFeature_SoupXRNA < rna.max.features &

            percent.mt < mt.threshold[[sample$name]]

    )

    DefaultAssay(obj) <- "ATAC"

    obj <- FindTopFeatures(

        obj,

        min.cutoff = lsi.min.cutoff

    )

    anchor.features <- VariableFeatures(obj)

    obj <- RunTFIDF(obj)

    obj <- RunSVD(obj)

    sample$object <- obj

    message(

        sample$name,

        " : ",

        ncol(obj),

        " cells retained"

    )

    return(sample)

}


detect_doublets <- function(sample){

    message("Detecting doublets: ", sample$name)

    DefaultAssay(sample$object) <- "SoupXRNA"

    sce <- as.SingleCellExperiment(sample$object)

    sce <- scDblFinder::scDblFinder(sce)

    sample$object$DoubletClass <- sce$scDblFinder.class

    sample$object$DoubletScore <- sce$scDblFinder.score

    message(
        sample$name,
        " : ",
        sum(sample$object$DoubletClass == "doublet"),
        " doublets detected"
    )

    sample

}


remove_doublets <- function(sample){

    message("Removing doublets: ", sample$name)

    sample$object <- subset(

        sample$object,

        subset = DoubletClass == "singlet"

    )

    message(

        sample$name,

        " : ",

        ncol(sample$object),

        " singlets retained"

    )

    sample

}




rename_cells <- function(sample){

    message("Renaming cells: ", sample$name)

    obj <- sample$object

    condition <- sample.condition[[sample$name]]

    if(is.null(condition))
        stop("Unknown condition for ", sample$name)

    obj <- RenameCells(

        obj,

        new.names = paste0(
            condition,
            "_",
            sample$name,
            "_",
            colnames(obj)
        )

    )

    obj$sample <- sample$name

    obj$condition <- condition

    sample$object <- obj

    return(sample)

}




merge_samples <- function(...){

    samples <- list(...)

    objects <- lapply(
        samples,
        function(x) x$object
    )

    names(objects) <- vapply(
        samples,
        function(x) x$name,
        character(1)
    )

    message(
        "Merging ",
        length(objects),
        " samples (",
        paste(names(objects), collapse = ", "),
        ")"
    )

    obj <- objects[[1]]

    if(length(objects) > 1){

        obj <- merge(

            x = objects[[1]],

            y = objects[-1]

        )

    }

    message(
        ncol(obj),
        " cells after merge"
    )

    return(obj)

}




run_atac_integration <- function(
    merged.obj,
    sample.list
){

    message("Running ATAC integration")

    ############################################################
    ## LSI INDIVIDUEL DE CHAQUE SAMPLE
    ############################################################

    sample.list <- lapply(
        sample.list,
        function(x){

            DefaultAssay(x) <- "ATAC"

            x <- FindTopFeatures(
                x,
                min.cutoff = 20
            )

            x <- RunTFIDF(
                x
            )

            x <- RunSVD(
                x
            )

            x
        }
    )

    ############################################################
    ## INTEGRATION RLSI
    ############################################################

    oplan <- future::plan()
    future::plan("sequential")

    integration.anchors <- FindIntegrationAnchors(

        object.list = sample.list,

        anchor.features = rownames(sample.list[[1]]),

        reduction = "rlsi",

        dims = 2:30

    )

    message(
        "Number of anchors found: ",
        nrow(integration.anchors@anchors)
    )

    ############################################################
    ## INTEGRATED LSI
    ############################################################

    DefaultAssay(merged.obj) <- "ATAC"

    merged.obj <- RunTFIDF(
        merged.obj
    )

    merged.obj <- FindTopFeatures(
        merged.obj,
        min.cutoff = 20
    )

    merged.obj <- RunSVD(
        merged.obj
    )

    lsi_for_harmony <- merged.obj[["lsi"]]

    merged.obj <- IntegrateEmbeddings(

        anchorset = integration.anchors,

        reductions = merged.obj[["lsi"]],

        new.reduction.name = "integrated_lsi",

        dims.to.integrate = 2:30,

        k.weight = 5

    )

    ############################################################
    ## HARMONY SUR LSI NON INTEGRE
    ############################################################

    merged.obj[["lsi"]] <- lsi_for_harmony

    merged.obj$batch <- ifelse(
        merged.obj$sample %in% c(
            "PG2",
            "PG6",
            "PG13",
            "PG24"
        ),
        "old",
        "new"
    )

    merged.obj <- RunHarmony(
        object = merged.obj,
        group.by.vars = "batch",
        reduction = "lsi",
        reduction.save = "harmony_lsi",
        project.dim = FALSE,
        verbose = TRUE
    )

    ############################################################
    ## UMAP SUR HARMONY LSI
    ############################################################

    future::plan(oplan)

    merged.obj <- RunUMAP(

        object = merged.obj,

        reduction = "harmony_lsi",

        dims = 2:30,

        reduction.name = "integrated_atac_umap"

    )

    return(merged.obj)

}






plot_atac_integration <- function(obj){

    out.dir <- file.path(qc.dir, "integration")

    dir.create(
        out.dir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    p <- DimPlot(

        obj,

        reduction = "integrated_atac_umap",

        group.by = "condition"

    ) +
    ggtitle("Integrated ATAC")

    ggsave(

        filename = file.path(
            out.dir,
            "ATAC_integration.pdf"
        ),

        plot = p,

        width = 6,
        height = 6

    )

    print(p)

}




call_peaks <- function(
    obj,
    annotation
){

    message("Calling peaks with MACS3")

    DefaultAssay(obj) <- "ATAC"

    peaks <- CallPeaks(

        object = obj,

        macs2.path = macs2.path

    )

    message(length(peaks), " peaks called")

    counts <- FeatureMatrix(

        fragments = Fragments(obj),

        features = peaks,

        cells = colnames(obj)

    )

    obj[["peaks"]] <- CreateChromatinAssay(

        counts = counts,

        fragments = Fragments(obj),

        annotation = annotation

    )

    message(

        nrow(obj[["peaks"]]),
        " quantified peaks"

    )

    return(obj)

}


process_peaks <- function(obj){

    message("Processing peak assay")

    DefaultAssay(obj) <- "peaks"

    obj <- FindTopFeatures(

        obj,

        min.cutoff = 5

    )

    obj <- RunTFIDF(obj)

    obj <- RunSVD(obj)

    obj <- RunUMAP(

        object = obj,

        reduction = "lsi",

        dims = 2:30,  #2:40 avant

        reduction.name = "atac_umap"

    )

    return(obj)

}



prepare_rna <- function(obj){

    message("Preparing RNA assay")

    DefaultAssay(obj) <- "SoupXRNA"

    mito.features <- c(
        "ND1",
        "ND2",
        "COX1",
        "COX2",
        "ATP8",
        "ATP6",
        "COX3",
        "ND3",
        "ND4L",
        "ND4",
        "ND5",
        "ND6",
        "CYTB"
    )

message("[1/7] Calculating mitochondrial percentage")

    obj$percent.mt <- PercentageFeatureSet(
        obj,
        features = mito.features
    )

    return(obj)

}




plot_rna_qc <- function(obj){

    message("RNA QC")

    out.dir <- file.path(
        qc.dir,
        "rna_qc"
    )

    dir.create(
        out.dir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    keep <- with(
        obj@meta.data,
        nFeature_SoupXRNA > 200 &
        nFeature_SoupXRNA < 4500 &
        percent.mt < 3 #3
    )

    message(
        "Cells before filtering: ",
        ncol(obj)
    )

    message(
        "Cells after historical filter: ",
        sum(keep)
    )

    message(
        "Cells removed: ",
        sum(!keep)
    )

    message("Per condition:")

    print(table(obj$condition))

    print(table(obj$condition[keep]))

    pdf(
        file.path(
            out.dir,
            "RNA_QC.pdf"
        ),
        width = 12,
        height = 10
    )

    print(
        VlnPlot(
            obj,
            features = c(
                "nFeature_SoupXRNA",
                "nCount_SoupXRNA",
                "percent.mt"
            ),
            ncol = 3,
            pt.size = 0
        )
    )

    print(
        FeaturePlot(
            obj,
            reduction = "integrated_atac_umap",
            features = "percent.mt",
            order = TRUE
        )
    )

    dev.off()

    invisible(keep)

}



filter_rna <- function(obj){

    message("Filtering RNA")

    keep <- with(
        obj@meta.data,
        nFeature_SoupXRNA > 200 &
        nFeature_SoupXRNA < 4500 &
        percent.mt < 3 #3
    )

    before <- ncol(obj)

    obj <- subset(
        obj,
        cells = colnames(obj)[keep]
    )

    message(
        before,
        " -> ",
        ncol(obj),
        " cells retained"
    )

    return(obj)

}


###############################################################################
## MAIN
###############################################################################

library(EnsDb.Rnorvegicus.v110)
library(BSgenome.Rnorvegicus.NCBI.rn7.2)

annotation <- GetGRangesFromEnsDb(
    ensdb = EnsDb.Rnorvegicus.v110
)

#GenomeInfoDb::seqlevelsStyle(annotation) <- "UCSC"

seqlevels(annotation) <- c(
    paste0("chr", 1:20),
    "chrX",
    "chrY",
    "chrM"
)


pg2  <- load_sample("PG2")
pg6  <- load_sample("PG6")
pg13 <- load_sample("PG13")
pg24 <- load_sample("PG24")
pg1  <- load_sample("PG1")
pg5  <- load_sample("PG5")
pg9  <- load_sample("PG9")
pg11  <- load_sample("PG11")
pg32  <- load_sample("PG32")
pg38  <- load_sample("PG38")
pgo1  <- load_sample("PGo1")
pgo6  <- load_sample("PGo6")





samples <- list(
    pg2,
    pg6,
    pg13,
    pg24,
    pg1,
    pg5,
    pg9,
    pg11,
    pg32,
    pg38,
    pgo1,
    pgo6
)


names(samples) <- c(
    "PG2",
    "PG6",
    "PG13",
    "PG24",
    "PG1",
    "PG5",
    "PG9",
    "PG11",
    "PG32",
    "PG38",
    "PGo1",
    "PGo6"
)

# ============================================================
# MACS3 PEAK CALLING — 12 SAMPLES
# ============================================================

sample_peaks <- lapply(
    names(cellranger.dir),
    function(s) {
        call_macs3_sample(
            sample = s,
            sample_dir = cellranger.dir[[s]]
        )
    }
)

names(sample_peaks) <- names(cellranger.dir)

peak_numbers <- lengths(sample_peaks)

print(
    data.frame(
        sample = names(sample_peaks),
        n_peaks = as.integer(peak_numbers)
    )
)

# ============================================================
# CONSENSUS PEAKS
# ============================================================

all_sample_peaks <- Reduce(
    c,
    sample_peaks
)

message(
    "Individual MACS3 peaks: ",
    length(all_sample_peaks)
)

consensus.peaks <- reduce(
    all_sample_peaks,
    ignore.strand = TRUE
)

consensus.peaks <- consensus.peaks[
    width(consensus.peaks) > 20 &
    width(consensus.peaks) < 10000
]

consensus.peaks <- keepStandardChromosomes(
    consensus.peaks,
    pruning.mode = "coarse"
)

consensus.peaks <- sort(consensus.peaks)

message(
    "Consensus MACS3 peaks: ",
    length(consensus.peaks)
)


for(i in seq_along(samples)){

    samples[[i]] <- run_soupx(samples[[i]])

    plot_rho_qc(samples[[i]])

    samples[[i]] <- apply_soupx(samples[[i]])

    samples[[i]] <- create_seurat(
        samples[[i]],
        consensus.peaks,
        annotation
    )

    plot_soupx_qc(samples[[i]])

    samples[[i]] <- plot_cell_qc(samples[[i]])

    samples[[i]] <- filter_cells(samples[[i]])

    samples[[i]] <- detect_doublets(samples[[i]])

    samples[[i]] <- remove_doublets(samples[[i]])

    samples[[i]] <- rename_cells(samples[[i]])

    gc()

}

pg2  <- samples[[1]]
pg6  <- samples[[2]]
pg13 <- samples[[3]]
pg24 <- samples[[4]]
pg1  <- samples[[5]]
pg5  <- samples[[6]]
pg9  <- samples[[7]]
pg11 <- samples[[8]]
pg32 <- samples[[9]]
pg38 <- samples[[10]]
pgo1 <- samples[[11]]
pgo6 <- samples[[12]]



obj <- merge_samples(
    pg2,
    pg6,
    pg13,
    pg24,
    pg1,
    pg5,
    pg9,
    pg11,
    pg32,
    pg38,
    pgo1,
    pgo6

)

obj <- run_atac_integration(

    merged.obj = obj,

    sample.list = list(
        pg2$object,
        pg6$object,
        pg13$object,
        pg24$object,
        pg1$object,
        pg5$object,
        pg9$object,
        pg11$object,
        pg32$object,
        pg38$object,
        pgo1$object,
        pgo6$object
    )

)

plot_atac_integration(obj)


obj <- RegionStats(
    object = obj,
    genome = BSgenome.Rnorvegicus.NCBI.rn7.2
)

obj <- prepare_rna(obj)

plot_rna_qc(obj)

obj <- filter_rna(obj)

saveRDS(
    obj,
    file.path(
        data.dir,
        "PGintegrated73_preWNN-harmony.rds"
    )
)

message("Pre-WNN object saved.")


process_rna <- function(obj){

message("[3/7] Normalizing RNA")

    DefaultAssay(obj) <- "SoupXRNA"

    obj <- NormalizeData(obj)

message("[4/7] Running SCTransform (long step)")

    obj <- SCTransform(
        obj,
        vars.to.regress = "percent.mt",
        verbose = FALSE
    )

    DefaultAssay(obj) <- "SoupXRNA"

message("[5/7] Finding variable features")

    obj <- FindVariableFeatures(obj)

message("[6/7] Running PCA")

oplan <- future::plan()
future::plan("sequential")

obj <- ScaleData(
        obj,
        verbose = FALSE
    )

    obj <- RunPCA(obj)


############################################################
## OPTIONAL : Harmony on RNA PCA
############################################################

 obj$batch <- ifelse(
     obj$sample %in% c("PG2","PG6","PG13","PG24"),
     "old",
     "new"
 )

 obj <- RunHarmony(
     object = obj,
     group.by.vars = "batch",
     reduction = "pca",
     reduction.save = "harmony_pca"
 )

############################################################



message("[7/7] Computing RNA UMAP")

    out.dir <- file.path(qc.dir, "rna")
    dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

    pdf(
        file.path(out.dir, "RNA_elbow.pdf"),
        width = 5,
        height = 5
    )
    print(ElbowPlot(obj))
    dev.off()

    obj <- RunUMAP(
        object = obj,
        reduction = "pca",
        assay = "SCT",
        reduction.name = "RNA_umap",
        dims = 2:30,
        verbose = TRUE
    )

    pdf(
        file.path(out.dir, "RNA_umap_condition.pdf"),
        width = 6,
        height = 6
    )
    print(
        DimPlot(
            obj,
            reduction = "RNA_umap",
            group.by = "condition"
        )
    )
    dev.off()

    pdf(
        file.path(out.dir, "RNA_umap_sample.pdf"),
        width = 6,
        height = 6
    )
    print(
        DimPlot(
            obj,
            reduction = "RNA_umap",
            group.by = "sample"
        )
    )
    dev.off()

    obj$sample_Id <- obj$orig.ident
    obj$orig.ident <- obj$condition

    return(obj)

}



process_multimodal <- function(obj){

    message("Building multimodal graph")


cat("PCA:", ncol(Embeddings(obj, "pca")), "\n")
cat("integrated_lsi:", ncol(Embeddings(obj, "integrated_lsi")), "\n")
cat("dims RNA:", paste(1:30, collapse = ","), "\n")
cat("dims ATAC:", paste(2:40, collapse = ","), "\n")


    obj <- FindMultiModalNeighbors(

        object = obj,

        reduction.list = list(
            "harmony_pca", #pca   #harmony_pca
            "harmony_lsi"  #harmony_lsi  #integrated_lsi  #lsi
        ),

        dims.list = list(
            1:30,
            2:30
        ),

        modality.weight.name = c(
    "RNA.weight",
    "peaks.weight"
),

        verbose = TRUE

    )

    obj <- RunUMAP(

        object = obj,

        nn.name = "weighted.nn",

        reduction.name = "multimodal_umap",

        assay = "RNA",

        return.neighbor = TRUE,

        verbose = TRUE

    )

    out.dir <- file.path(
        qc.dir,
        "multimodal"
    )

    dir.create(
        out.dir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    pdf(
        file.path(
            out.dir,
            "multimodal_umap_condition_perc-mt3_harm-lsi2-30_harmPCA1-30.pdf"
        ),
        width = 6,
        height = 6
    )

    print(
        DimPlot(
            obj,
            reduction = "multimodal_umap",
            group.by = "condition"
        )
    )

    dev.off()

    pdf(
        file.path(
            out.dir,
            "multimodal_umap_sample_perc-mt3_harm-lsi2-30_harmPCA1-30.pdf"
        ),
        width = 6,
        height = 6
    )

    print(
        DimPlot(
            obj,
            reduction = "multimodal_umap",
            group.by = "sample"
        )
    )

    dev.off()

    return(obj)

}


plot_integrated_qc <- function(obj){

    message("Integrated QC")

    DefaultAssay(obj) <- "SoupXRNA"

    out.dir <- file.path(
        qc.dir,
        "integrated_qc"
    )

    dir.create(
        out.dir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    pdf(
        file.path(
            out.dir,
            "Integrated_QC.pdf"
        ),
        width = 16,
        height = 14
    )

    print(
        FeaturePlot(
            obj,
            reduction = "multimodal_umap",
            features = c(
                "percent.mt",
                "nCount_SoupXRNA",
                "nFeature_SoupXRNA",
                "nCount_ATAC",
                "TSS.enrichment",
                "nucleosome_signal"
            ),
            order = TRUE,
            ncol = 3
        )
    )

    print(
        FeatureScatter(
            obj,
            feature1 = "nCount_SoupXRNA",
            feature2 = "nCount_ATAC"
        )
    )

    print(
        FeatureScatter(
            obj,
            feature1 = "nFeature_SoupXRNA",
            feature2 = "nCount_ATAC"
        )
    )

    print(
        FeatureScatter(
            obj,
            feature1 = "nCount_SoupXRNA",
            feature2 = "percent.mt"
        )
    )

    print(
        FeaturePlot(
            obj,
            reduction = "integrated_atac_umap",
            features = soupx.markers,
            order = TRUE,
            ncol = 2
        )
    )

    dev.off()

    invisible(obj)

}




cluster_multimodal <- function(
    obj,
    resolution = 0.4
){

    message("Clustering multimodal graph")

    out.dir <- file.path(
        qc.dir,
        "clusters"
    )

    dir.create(
        out.dir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    obj <- FindClusters(

        object = obj,

        graph = "wknn",

        resolution = resolution,

        verbose = FALSE,

        random.seed = 1234

    )

obj$cluster_id <- as.character(obj$seurat_clusters)

    p <- DimPlot(

        obj,

        reduction = "multimodal_umap",

        group.by = "seurat_clusters",

        label = TRUE

    ) +
    ggtitle(
        paste0(
            "Multimodal clusters (resolution = ",
            resolution,
            ")"
        )
    )

    ggsave(

        filename = file.path(
            out.dir,
            "Multimodal_clusters.pdf"
        ),

        plot = p,

        width = 7,

        height = 6

    )

    print(p)

    DefaultAssay(obj) <- "SoupXRNA"

    marker.genes <- c(

    "Gh1",
    "Prl",
    "Lhb",
    "Fshb",
    "Cga",
    "Tshb",
    "Pomc",
    "S100b",
    "Prop1",
    "Plvap",
    "Emcn",
    "Tyrobp",
    "Arhgap15",
    "Pde5a",
    "Pdgfrb",
    "Adamts9",
    "Rfx4",
    "Neurod1",
    "Pax7",
    "Pou1f1",
    "S100g",
    "Ghrhr",
    "Pappa2",
    "Dlk1",
    "Agtr1b",
    "Ralyl",
    "Gria4"

)

pdf(

    file.path(
        out.dir,
        "Cluster_marker_featureplots_perc-mt3_harm-lsi2-30_harmPCA1-30.pdf"
    ),

    width = 10,
    height = 10

)

for(i in seq(1, length(marker.genes), by = 4)){

    genes <- marker.genes[
        i:min(i + 3, length(marker.genes))
    ]

    print(

        FeaturePlot(

            obj,

            reduction = "multimodal_umap",

            features = genes,

            order = TRUE,

            max.cutoff = "q90",

            ncol = 2

        )

    )

}

dev.off()

    return(obj)

}

annotate_clusters <- function(obj){

    message("Annotating clusters")

    obj$cluster_id <- unname(
        cluster.annotation[
            as.character(obj$seurat_clusters)
        ]
    )

    p <- DimPlot(

        obj,

        reduction = "multimodal_umap",

        group.by = "cluster_id",

        label = TRUE,

        repel = TRUE

    ) +
    ggtitle("Annotated multimodal clusters")

    ggsave(

        filename = file.path(
            qc.dir,
            "clusters",
            "Annotated_clusters.pdf"
        ),

        plot = p,

        width = 8,

        height = 7

    )

    return(obj)

}





###############################################################################
## MAIN
###############################################################################

obj <- readRDS(
    file.path(
        data.dir,
        "PGintegrated73_preWNN-harmony.rds"
    )
)


gc()

message(
    format(
        object.size(obj),
        units = "GB"
    )
)

obj <- process_rna(obj)

obj <- process_multimodal(obj)

obj <- cluster_multimodal(
    obj,
    resolution = 0.4
)

plot_integrated_qc(obj)

obj <- annotate_clusters(obj)

saveRDS(

    obj,

    file.path(
        data.dir,
        "PGintegrated73.rebuilt-harm.rds"
    )

)

message("Pipeline completed successfully.")

