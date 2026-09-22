
suppressPackageStartupMessages({

  library(Seurat)
  library(Signac)
  library(scMEGA)
  library(mgcv)
  library(SummarizedExperiment)
  library(Matrix)

  library(data.table)
  library(dplyr)

})

options(stringsAsFactors = FALSE)


trajectory_prefix <- "pt_slingshot2d_"

cells_per_bin <- 10

gam_qvalue <- 0.01

k_mad_filter <- 0.75

k_mad_filter_atac <- 0.2

gam_k <- 5

velocity_r2_min <- 0.20

atac_smoothing_window <- 5

atac_r2_min <- 0.15


dir.create(
  "results/pseudotime_DE",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/pseudotime_DE/RNA",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/pseudotime_DE/Peaks",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/pseudotime_DE/Velocity",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/pseudotime_DE/objects",
  recursive = TRUE,
  showWarnings = FALSE
)


cat(
  "\n====================\n",
  "Chargement objet\n",
  "====================\n"
)

integrated73 <- readRDS(
  "results/objects/PGintegrated73.reclustered.linkpeaks.trajectory.rds"
)

cat(
  "Cellules :",
  ncol(integrated73),
  "\n"
)

cat(
  "Genes :",
  nrow(integrated73[["SoupXRNA"]]),
  "\n"
)


clusters <- "Cluster_G"

cat(
  "Clusters :",
  paste(clusters, collapse = ", "),
  "\n"
)


Dynamic_results <- list()

DE_results <- list()

Velocity_results <- list()

GetTrajectory_adaptiveBins <- function(
    object,
    trajectory.name = "Trajectory",
    assay,
    slot = "data",
    cells_per_bin = 10,
    log2Norm = TRUE,
    scaleTo = 10000
){

  if (is.null(assay) || !assay %in% names(object@assays))
    stop("Please provide a valid assay")

  if (!trajectory.name %in% colnames(object@meta.data))
    stop("Trajectory not found in meta.data")

  trajectory <- object@meta.data[[trajectory.name]]

  if (!is.numeric(trajectory))
    stop("Trajectory must be numeric")

  keep <- is.finite(trajectory)
  trajectory <- trajectory[keep]
  cells <- colnames(object)[keep]

  n_cells <- length(cells)

  if (n_cells < cells_per_bin * 2)
    stop("Not enough cells to build trajectory bins")

  nbins <- floor(n_cells / cells_per_bin)

  if (nbins < 2)
    stop("nbins < 2")

  breaks <- seq(0, 100, length.out = nbins + 1)

  bin_id <- cut(
    trajectory,
    breaks = breaks,
    include.lowest = TRUE,
    labels = FALSE
  )

  names(bin_id) <- cells

  data.use <- LayerData(
    object,
    assay = assay,
    layer = slot
  )

  groupMat <- sapply(
    split(cells, bin_id),
    function(cell_names){
      Matrix::rowMeans(
        data.use[, cell_names, drop = FALSE]
      )
    }
  )

  colnames(groupMat) <- paste0(
    "Bin_",
    seq_len(ncol(groupMat))
  )

  if (!is.null(scaleTo)) {

    if (any(groupMat < 0)) {

      message(
        "Values < 0 detected, skipping depth normalization"
      )

    } else {

      groupMat <- t(
        t(groupMat) /
        colSums(groupMat)
      ) * scaleTo

    }
  }

  if (log2Norm) {

    if (any(groupMat < 0)) {

      message(
        "Values < 0 detected, skipping log2 normalization"
      )

    } else {

      groupMat <- log2(groupMat + 1)

    }
  }


seTrajectory <- SummarizedExperiment(

  assays = SimpleList(

    mat = as.matrix(
      groupMat
    )

  )

)

return(
  seTrajectory
)

}

cat(
  "\nPARTIE 1 CHARGÉE\n"
)


load_velocity_stochastic <- function(
    cluster
){

  velocity_path <- file.path(
    "results",
    "velocyto",
    "exports",
    cluster,
    "velocity_stochastic_all.csv"
  )

  if(
    !file.exists(velocity_path)
  ){

    stop(
      "Fichier velocity absent : ",
      velocity_path
    )

  }

  cat(
    "Lecture velocity stochastic :",
    velocity_path,
    "\n"
  )

  velocity <- fread(
    velocity_path,
    data.table = FALSE
  )


  barcode <- velocity[[1]]

  velocity <- velocity[
    ,
    -1,
    drop = FALSE
  ]

  rownames(
    velocity
  ) <- barcode



  velocity <- as.matrix(
    velocity
  )

  mode(
    velocity
  ) <- "numeric"


  cat(
    "Cellules velocity :",
    nrow(velocity),
    "\n"
  )

  cat(
    "Gènes velocity :",
    ncol(velocity),
    "\n"
  )

  cat(
    "Première cellule :",
    rownames(velocity)[1],
    "\n"
  )

  cat(
    "Premier gène :",
    colnames(velocity)[1],
    "\n"
  )

  velocity[
    !is.finite(velocity)
  ] <- NA_real_

  return(
    velocity
  )

}


cat(
  "\n====================\n",
  "Construction trajectoires\n",
  "====================\n"
)

for(cl in clusters){

  cat(
    "\n--------------------\n"
  )

  cat(
    "Cluster:",
    cl,
    "\n"
  )


  traj_col <- paste0(
    trajectory_prefix,
    cl
  )

  if(
    !(traj_col %in%
        colnames(
          integrated73@meta.data
        ))
  ){

    cat(
      "Trajectoire absente:",
      traj_col,
      "\n"
    )

    next

  }

  cells <- colnames(
    integrated73
  )[
    integrated73$traj_cluster == cl
  ]

  if(
    length(cells) < 50
  ){

    cat(
      "skip (<50 cellules)\n"
    )

    next

  }

  traj_vec <- integrated73@meta.data[
    cells,
    traj_col
  ]

  names(traj_vec) <- cells

  keep <- is.finite(
    traj_vec
  )

  if(
    sum(keep) < 50
  ){

    cat(
      "skip (<50 cellules avec trajectoire)\n"
    )

    next

  }


  obj_tmp <- integrated73[
    ,
    names(traj_vec)[keep]
  ]

  traj_vec <- traj_vec[keep]

  traj_vec <- traj_vec -
    min(
      traj_vec,
      na.rm = TRUE
    )

  traj_vec <- traj_vec /
    max(
      traj_vec,
      na.rm = TRUE
    ) * 100

  obj_tmp$Trajectory_tmp <- traj_vec

  cat(
    "Cellules retenues:",
    length(traj_vec),
    "\n"
  )

DefaultAssay(obj_tmp) <- "SoupXRNA"

obj_rna <- obj_tmp

obj_rna <- NormalizeData(
  obj_rna,
  verbose = FALSE
)


trajRNA <- GetTrajectory_adaptiveBins(

  object = obj_rna,

  trajectory.name =
    "Trajectory_tmp",

  assay = "SoupXRNA",

  slot = "data",

  cells_per_bin =
    cells_per_bin,

  log2Norm = TRUE

)

expr_mat <- assay(
  trajRNA,
  "mat"
)

colnames(expr_mat) <- paste0(
  "bin_",
  seq_len(
    ncol(expr_mat)
  )
)

cat(
  "RNA :",
  nrow(expr_mat),
  "genes x",
  ncol(expr_mat),
  "bins\n"
)


DefaultAssay(obj_tmp) <- "peaks"

obj_peak <- obj_tmp

obj_peak <- RunTFIDF(
  obj_peak,
  verbose = FALSE
)

obj_peak <- FindTopFeatures(
  obj_peak,
  min.cutoff = "q0"
)

obj_peak <- RunSVD(
  obj_peak,
  verbose = FALSE
)

trajATAC <- GetTrajectory_adaptiveBins(

  object = obj_peak,

  trajectory.name =
    "Trajectory_tmp",

  assay = "peaks",

  slot = "data",

  cells_per_bin =
    cells_per_bin,

  log2Norm = FALSE

)

peak_mat <- assay(
  trajATAC,
  "mat"
)

colnames(peak_mat) <- paste0(
  "bin_",
  seq_len(
    ncol(peak_mat)
  )
)

cat(
  "ATAC :",
  nrow(peak_mat),
  "peaks x",
  ncol(peak_mat),
  "bins\n"
)


  Dynamic_results[[cl]] <- list(

    expr_mat = expr_mat,

    peak_mat = peak_mat,

    trajectory = traj_vec,

    traj_col = traj_col

  )

}

cat(
  "\nPARTIE 2 TERMINÉE\n"
)

cat(
  "Clusters traités :\n"
)

print(
  names(
    Dynamic_results
  )
)


cat(
  "\n====================\n",
  "DEG RNA\n",
  "====================\n"
)

gene_pattern <- c(
  "^LOC",
  "^NEWGENE",
  "^RGD",
  "^ENSRN",
  "^AC[0-9]"
)

for(cl in names(Dynamic_results)){

  cat(
    "\n--------------------\n"
  )

  cat(
    "RNA:",
    cl,
    "\n"
  )

  expr_mat <- Dynamic_results[[cl]]$expr_mat


  expr_mat <- expr_mat[

    !grepl(

      paste(
        gene_pattern,
        collapse = "|"
      ),

      rownames(expr_mat)

    ),

    ,
    drop = FALSE

  ]

  cat(
    "Genes après filtre annotation:",
    nrow(expr_mat),
    "\n"
  )


sdv <- apply(
  expr_mat,
  1,
  sd,
  na.rm = TRUE
)

thr <- median(
  sdv,
  na.rm = TRUE
) +
  k_mad_filter *
  mad(
    sdv,
    na.rm = TRUE
  )

expr_filt <- expr_mat[
  sdv > thr,
  ,
  drop = FALSE
]

cat(
  "Seuil MAD RNA:",
  thr,
  "\n"
)

cat(
  "Genes après filtre MAD:",
  nrow(expr_filt),
  "\n"
)

##########################################################
## GAM
##########################################################

cat(
  "Genes testés par GAM:",
  nrow(expr_filt),
  "\n"
)



time_bins <- seq_len(
  ncol(expr_filt)
)

pvals <- apply(

  expr_filt,

  1,

  function(e){

    df <- data.frame(

      expression = as.numeric(e),

      pseudotime = time_bins

    )

    out <- try(

      gam(
        expression ~
          s(
            pseudotime,
            k = gam_k,
            bs = "tp"
          ),
        data = df,
        method = "REML"
      ),

      silent = TRUE

    )

    if(
      inherits(
        out,
        "try-error"
      )
    ){

      return(1)

    }

    sm <- summary(
      out
    )

    ## p-value du terme spline

    p <- sm$s.table[
      1,
      "p-value"
    ]

    if(
      !is.finite(p)
    ){

      return(1)

    }

    p

  }

)

qvals <- p.adjust(
  pvals,
  method = "BH"
)

keep <- qvals <
  gam_qvalue

  expr_keep <- expr_filt[
    keep,
    ,
    drop = FALSE
  ]

  RNA_dyn <- rownames(
    expr_keep
  )

  ##########################################################
  ## TABLE
  ##########################################################

  RNA_table <- data.table(

    gene = RNA_dyn,

    amplitude = apply(

      expr_keep,

      1,

      function(x){

        quantile(
          x,
          0.95,
          na.rm = TRUE
        ) -

          quantile(
            x,
            0.05,
            na.rm = TRUE
          )

      }

    ),

    p_value = pvals[keep],

    q_value = qvals[keep]

  )

  setorder(
    RNA_table,
    q_value
  )

  cat(
    "RNA dyn:",
    nrow(RNA_table),
    "\n"
  )

  ##########################################################
  ## EXPORT CSV
  ##########################################################

  fwrite(

    RNA_table,

    file.path(

      "results/pseudotime_DE/RNA",

      paste0(
        cl,
        "_RNA.csv"
      )

    )

  )

  DE_results[[cl]] <- list(

    RNA_dyn = RNA_dyn,

    RNA_table = RNA_table,

    expr_keep = expr_keep

  )

}

cat(
  "\nPARTIE 3 TERMINÉE\n"
)

cat(
  "Clusters RNA disponibles :\n"
)

print(
  names(
    DE_results
  )
)

############################################################
## PART 4
## DAE ATAC
############################################################

cat(
  "\n====================\n",
  "DAR ATAC\n",
  "====================\n"
)

for(cl in names(Dynamic_results)){

  cat(
    "\n--------------------\n"
  )

  cat(
    "ATAC:",
    cl,
    "\n"
  )


  peak_mat <- Dynamic_results[[cl]]$peak_mat

  cat(
    "Peaks avant filtre MAD:",
    nrow(peak_mat),
    "\n"
  )

  cat(
    "Bins:",
    ncol(peak_mat),
    "\n"
  )


  peak_sdv <- apply(
    peak_mat,
    1,
    sd,
    na.rm = TRUE
  )

  peak_thr <- median(
    peak_sdv,
    na.rm = TRUE
  ) +
    k_mad_filter_atac *
    mad(
      peak_sdv,
      na.rm = TRUE
    )

  peak_filt <- peak_mat[
    peak_sdv > peak_thr,
    ,
    drop = FALSE
  ]

  cat(
    "Seuil MAD ATAC:",
    peak_thr,
    "\n"
  )

  cat(
    "k MAD:",
    k_mad_filter_atac,
    "\n"
  )

  cat(
    "Peaks après filtre MAD:",
    nrow(peak_filt),
    "\n"
  )

  cat(
    "\nSmoothing ATAC avec scMEGA\n"
  )

  cat(
    "Méthode : scMEGA:::centerRollMean\n"
  )

  cat(
    "Fenêtre :",
    atac_smoothing_window,
    "bins\n"
  )

  peak_smooth <- t(
    apply(
      peak_filt,
      1,
      function(x){

        scMEGA:::centerRollMean(
          v = x,
          k = atac_smoothing_window
        )

      }
    )
  )

  peak_smooth <- as.matrix(
    peak_smooth
  )

  rownames(peak_smooth) <-
    rownames(peak_filt)

  colnames(peak_smooth) <-
    colnames(peak_filt)

  cat(
    "Matrice ATAC après smoothing :",
    nrow(peak_smooth),
    "peaks x",
    ncol(peak_smooth),
    "bins\n"
  )


  cat(
    "\nGAM ATAC\n"
  )

  time_bins <- seq_len(
    ncol(peak_smooth)
  )


  pvals <- apply(

    peak_smooth,

    1,

    function(e){

      df <- data.frame(

        accessibility =
          as.numeric(e),

        pseudotime =
          time_bins

      )

      fit <- try(

        gam(
          accessibility ~
            s(
              pseudotime,
              k = gam_k,
              bs = "tp"
            ),
          data = df,
          method = "REML"
        ),

        silent = TRUE

      )

      if(
        inherits(
          fit,
          "try-error"
        )
      ){

        return(1)

      }

      sm <- summary(
        fit
      )

      p <- sm$s.table[
        1,
        "p-value"
      ]

      if(
        !is.finite(p)
      ){

        return(1)

      }

      p

    }

  )


  r2vals <- apply(

    peak_smooth,

    1,

    function(e){

      df <- data.frame(

        accessibility =
          as.numeric(e),

        pseudotime =
          time_bins

      )

      fit <- try(

        gam(
          accessibility ~
            s(
              pseudotime,
              k = gam_k,
              bs = "tp"
            ),
          data = df,
          method = "REML"
        ),

        silent = TRUE

      )

      if(
        inherits(
          fit,
          "try-error"
        )
      ){

        return(0)

      }

      r2 <- summary(
        fit
      )$r.sq

      if(
        !is.finite(r2)
      ){

        return(0)

      }

      r2

    }

  )


  qvals <- p.adjust(
    pvals,
    method = "BH"
  )


  cat(
    "Peaks testés par GAM:",
    length(pvals),
    "\n"
  )

  cat(
    "ATAC q < 0.05:",
    sum(
      qvals < gam_qvalue,
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "ATAC R² > 0.20:",
    sum(
      r2vals > atac_r2_min,
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "ATAC q < 0.05 & R² > 0.20:",
    sum(
      qvals < gam_qvalue &
      r2vals > atac_r2_min,
      na.rm = TRUE
    ),
    "\n"
  )


  keep <-
    qvals < gam_qvalue &
    r2vals > atac_r2_min

  peak_keep <- peak_smooth[
    keep,
    ,
    drop = FALSE
  ]

  Peak_dyn <- rownames(
    peak_keep
  )


  Peak_table <- data.table(

    peak = Peak_dyn,

    amplitude = apply(

      peak_keep,

      1,

      function(x){

        quantile(
          x,
          0.95,
          na.rm = TRUE
        ) -

          quantile(
            x,
            0.05,
            na.rm = TRUE
          )

      }

    ),

    p_value = pvals[
      keep
    ],

    q_value = qvals[
      keep
    ],

    r2 = r2vals[
      keep
    ]

  )

  setorder(
    Peak_table,
    q_value,
    -r2
  )

  cat(
    "\nATAC dyn :",
    nrow(Peak_table),
    "\n"
  )


  fwrite(

    Peak_table,

    file.path(

      "results/pseudotime_DE/Peaks",

      paste0(
        cl,
        "_Peaks.csv"
      )

    )

  )

  if(
    !cl %in% names(DE_results)
  ){

    DE_results[[cl]] <- list()

  }

  DE_results[[cl]]$Peak_dyn <-
    Peak_dyn

  DE_results[[cl]]$Peak_table <-
    Peak_table

  DE_results[[cl]]$peak_keep <-
    peak_keep

}

cat(
  "\nPARTIE 4 TERMINÉE\n"
)

cat(
  "Clusters ATAC disponibles :\n"
)

print(
  names(
    DE_results
  )
)



cat(
  "\n====================\n",
  "VELOCITY STOCHASTIC\n",
  "====================\n"
)

for(cl in clusters){

  cat(
    "\n--------------------\n"
  )

  cat(
    "Velocity:",
    cl,
    "\n"
  )


  cells <- colnames(
    integrated73
  )[
    integrated73$traj_cluster == cl
  ]


  traj_col <- paste0(
    trajectory_prefix,
    cl
  )

  if(
    !traj_col %in%
    colnames(
      integrated73@meta.data
    )
  ){

    stop(
      "Trajectoire absente : ",
      traj_col
    )

  }

  traj <- integrated73@meta.data[
    cells,
    traj_col
  ]

  names(traj) <- cells

  keep <- is.finite(
    traj
  )

  cells <- cells[
    keep
  ]

  traj <- traj[
    keep
  ]

  cat(
    "Cellules avec pseudotemps :",
    length(cells),
    "\n"
  )


  velocity <- load_velocity_stochastic(
    cl
  )


  common <- intersect(
    cells,
    rownames(velocity)
  )

  cat(
    "Cellules communes :",
    length(common),
    "\n"
  )

  if(
    length(common) < 50
  ){

    stop(
      "Moins de 50 cellules communes entre Seurat et velocity"
    )

  }

  velocity <- velocity[
    common,
    ,
    drop = FALSE
  ]

  traj <- traj[
    common
  ]

  ord <- order(
    traj
  )

  common <- common[
    ord
  ]

  traj <- traj[
    ord
  ]

  velocity <- velocity[
    common,
    ,
    drop = FALSE
  ]
velocity_cell <- velocity

  bins <- split(

    seq_along(traj),

    ceiling(
      seq_along(traj) /
      cells_per_bin
    )

  )

  velocity_mat <- sapply(

    bins,

    function(idx){

      colMeans(
        velocity[
          idx,
          ,
          drop = FALSE
        ],
        na.rm = TRUE
      )

    }

  )

  velocity_mat <- as.matrix(
    velocity_mat
  )


  colnames(
    velocity_mat
  ) <- paste0(
    "bin_",
    seq_len(
      ncol(velocity_mat)
    )
  )

  cat(
    "Velocity :",
    nrow(velocity_mat),
    "genes x",
    ncol(velocity_mat),
    "bins\n"
  )


velocity_annotation_keep <- !grepl(

  paste(
    gene_pattern,
    collapse = "|"
  ),

  rownames(velocity_mat)

)

cat(
  "Gènes Velocity avant filtre annotation:",
  nrow(velocity_mat),
  "\n"
)

velocity_mat <- velocity_mat[
  velocity_annotation_keep,
  ,
  drop = FALSE
]

cat(
  "Gènes Velocity après filtre annotation:",
  nrow(velocity_mat),
  "\n"
)


velocity_sdv <- apply(
  velocity_mat,
  1,
  sd,
  na.rm = TRUE
)

velocity_thr <- median(
  velocity_sdv,
  na.rm = TRUE
) +
  k_mad_filter *
  mad(
    velocity_sdv,
    na.rm = TRUE
  )

velocity_filt <- velocity_mat[
  velocity_sdv > velocity_thr,
  ,
  drop = FALSE
]

cat(
  "Seuil MAD Velocity:",
  velocity_thr,
  "\n"
)

cat(
  "Gènes velocity après filtre MAD:",
  nrow(velocity_filt),
  "\n"
)

cat(
  "Gènes velocity testés par GAM:",
  nrow(velocity_filt),
  "\n"
)


time_bins <- seq_len(
  ncol(velocity_filt)
)



pvals <- apply(

  velocity_filt,

  1,

  function(e){

    df <- data.frame(

      velocity = as.numeric(e),

      pseudotime = time_bins

    )

    fit <- try(

      gam(
        velocity ~
          s(
            pseudotime,
            k = gam_k,
            bs = "tp"
          ),
        data = df,
        method = "REML"
      ),

      silent = TRUE

    )

    if(
      inherits(
        fit,
        "try-error"
      )
    ){

      return(1)

    }

    sm <- summary(
      fit
    )

    p <- sm$s.table[
      1,
      "p-value"
    ]

    if(
      !is.finite(p)
    ){

      return(1)

    }

    p

  }

)


r2vals <- apply(

  velocity_filt,

  1,

  function(e){

    df <- data.frame(

      velocity = as.numeric(e),

      pseudotime = time_bins

    )

    fit <- try(

      gam(
        velocity ~
          s(
            pseudotime,
            k = gam_k,
            bs = "tp"
          ),
        data = df,
        method = "REML"
      ),

      silent = TRUE

    )

    if(
      inherits(
        fit,
        "try-error"
      )
    ){

      return(0)

    }

    r2 <- summary(
      fit
    )$r.sq

    if(
      !is.finite(r2)
    ){

      return(0)

    }

    r2

  }

)



qvals <- p.adjust(
  pvals,
  method = "BH"
)


cat(
  "Velocity q < 0.05 :",
  sum(
    qvals < gam_qvalue,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Velocity R² > 0.20 :",
  sum(
    r2vals > velocity_r2_min,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Velocity q < 0.05 & R² > 0.20 :",
  sum(
    qvals < gam_qvalue &
    r2vals > velocity_r2_min,
    na.rm = TRUE
  ),
  "\n"
)

keep <- (

  qvals < gam_qvalue

) & (

  r2vals > velocity_r2_min

)

velocity_keep <- velocity_filt[
  keep,
  ,
  drop = FALSE
]

Velocity_dyn <- rownames(
  velocity_keep
)



Velocity_table <- data.table(

  gene = Velocity_dyn,

  p_value = pvals[
    Velocity_dyn
  ],

  q_value = qvals[
    Velocity_dyn
  ],

  r2 = r2vals[
    Velocity_dyn
  ]

)

setorder(
  Velocity_table,
  q_value
)
cat(
  "Velocity dyn (q < ",
  gam_qvalue,
  " & R² > ",
  velocity_r2_min,
  ") : ",
  nrow(Velocity_table),
  "\n",
  sep = ""
)



  fwrite(

    Velocity_table,

    file.path(

      "results/pseudotime_DE/Velocity",

      paste0(
        cl,
        "_Velocity_stochastic.csv"
      )

    )

  )


  DE_results[[cl]]$Velocity_dyn <-
    Velocity_dyn

  DE_results[[cl]]$Velocity_table <-
    Velocity_table


Velocity_results[[cl]] <- list(

  velocity_cell =
    velocity_cell,

  velocity_mat =
    velocity_mat,

  velocity_filt =
    velocity_filt,

  p_value =
    pvals,

  q_value =
    qvals,

  Velocity_dyn =
    Velocity_dyn,

  Velocity_table =
    Velocity_table,

  trajectory =
    traj,

  traj_col =
    traj_col

)
}

cat(
  "\nPARTIE 5 TERMINÉE\n"
)

cat(
  "Clusters Velocity disponibles :\n"
)

print(
  names(
    Velocity_results
  )
)


cat(
  "\n====================\n",
  "Sauvegardes finales\n",
  "====================\n"
)


summary_list <- list()

for(cl in names(DE_results)){

  n_rna <- 0
  n_peak <- 0
  n_velocity <- 0

  if(
    "RNA_table" %in%
    names(DE_results[[cl]])
  ){

    n_rna <- nrow(
      DE_results[[cl]]$RNA_table
    )

  }

  if(
    "Peak_table" %in%
    names(DE_results[[cl]])
  ){

    n_peak <- nrow(
      DE_results[[cl]]$Peak_table
    )

  }

if(
  "Velocity_table" %in%
  names(DE_results[[cl]])
){

  n_velocity <- nrow(
    DE_results[[cl]]$Velocity_table
  )

}

  summary_list[[cl]] <- data.frame(

    cluster = cl,

    n_DEG = n_rna,

    n_DAR = n_peak,

    n_velocity = n_velocity,

    stringsAsFactors = FALSE

  )

}

DE_summary <- do.call(
  rbind,
  summary_list
)

rownames(
  DE_summary
) <- NULL



fwrite(

  DE_summary,

  file =
    "results/pseudotime_DE/DE_summary.csv"

)


saveRDS(

  Dynamic_results,

  file =
    "results/pseudotime_DE/objects/Dynamic_results.rds"

)


saveRDS(

  Velocity_results,

  file =
    "results/pseudotime_DE/objects/Velocity_results.rds"

)


saveRDS(

  DE_results,

  file =
    "results/pseudotime_DE/objects/DE_results.rds"

)



cat(
  "\nRésumé :\n"
)

print(
  DE_summary
)

cat(
  "\nClusters présents :\n"
)

print(
  names(DE_results)
)

cat(
  "\nObjets sauvegardés :\n"
)

cat(
  "results/pseudotime_DE/objects/Dynamic_results.rds\n"
)

cat(
  "results/pseudotime_DE/objects/DE_results.rds\n"
)

cat(
  "results/pseudotime_DE/DE_summary.csv\n"
)


cat(
  "\n========================================\n"
)

cat(
  "08_pseudotime_DEG_DAR.R terminé\n"
)

cat(
  "Trajectoire utilisée : ",
  trajectory_prefix,
  "\n",
  sep = ""
)

cat(
  "Clusters traités : ",
  length(names(DE_results)),
  "\n",
  sep = ""
)

cat(
  "========================================\n"
)