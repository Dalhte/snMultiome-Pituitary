
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(scMEGA)
library(slingshot)
library(SingleCellExperiment)
  library(dbscan)

  library(dplyr)
  library(tidyr)
  library(readr)
  library(mgcv)

  library(ggplot2)
})

set.seed(1234)


dir.create(
  "results/trajectory",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/trajectory/DBSCAN_outliers",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/trajectory/Slingshot2D",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "results/trajectory/Frises",
  recursive = TRUE,
  showWarnings = FALSE
)


integrated73 <- readRDS(
  "results/objects/PGintegrated73.reclustered.linkpeaks.rds"
)

DefaultAssay(integrated73) <- "SoupXRNA"

cat(
  "Objet chargé :",
  ncol(integrated73),
  "cellules\n"
)



red_2d <- "multimodal_umap_2D"

eps_2d <- 0.35
minPts <- 45
keep_n_cores <- 2

eta_group <- "orig.ident"

velocyto_root <- "results/velocyto/exports"

group_colors <- c(
  D2 = "#0000FF",
  PM = "#00FF00",
  PS = "#FF0000",
  E  = "#FFFF00"
)


cat("\nConstruction de traj_cluster...\n")

integrated73$traj_cluster <-
  as.character(integrated73$new_cluster)

############################################################
## Fusion L1/L2
############################################################

integrated73$traj_cluster[
  integrated73$traj_cluster %in%
    c("Cluster_L1", "Cluster_L2")
] <- "Cluster_L"

############################################################
## Fusion G1/G2
############################################################

integrated73$traj_cluster[
  integrated73$traj_cluster %in%
    c("Cluster_G1", "Cluster_G2")
] <- "Cluster_G"


clusters <- unique(
  integrated73$traj_cluster
)

clusters <- clusters[
  !is.na(clusters)
]

clusters <- sort(clusters)

cat("\nClusters trajectoire :\n")

print(clusters)

cat("\nEffectifs :\n")

print(
  table(integrated73$traj_cluster)
)


rescale01 <- function(x){

  rng <- range(
    x,
    na.rm = TRUE
  )

  if(diff(rng) == 0){

    return(
      rep(
        0.5,
        length(x)
      )
    )

  }

  (x - rng[1]) / diff(rng)
}


eta2_from_aov <- function(x, grp){

  df <- data.frame(
    x   = as.numeric(x),
    grp = as.factor(grp)
  )

  df <- df[
    complete.cases(df),
  ]

  if(nrow(df) < 30)
    return(NA_real_)

  if(nlevels(df$grp) < 2)
    return(NA_real_)

  fit <- aov(
    x ~ grp,
    data = df
  )

  tab <- summary(fit)[[1]]

  ss_eff <- tab[
    1,
    "Sum Sq"
  ]

  ss_tot <- sum(
    tab[, "Sum Sq"]
  )

  as.numeric(
    ss_eff / ss_tot
  )
}


compute_unified_time <- function(
    pt,
    latent
){

  df <- data.frame(
    pt      = as.numeric(pt),
    latent  = as.numeric(latent)
  )

  keep <- complete.cases(df)

  if(sum(keep) < 50){

    return(
      rep(
        NA_real_,
        length(pt)
      )
    )

  }

  theta <- rescale01(
    df$pt[keep]
  ) * 2 * pi

  gam_mod <- gam(

    df$latent[keep] ~

      s(
        theta,
        bs = "cc",
        k = 10
      ),

    knots = list(
      theta = c(
        0,
        2 * pi
      )
    ),

    method = "REML"
  )

  theta_all <- rescale01(pt) * 2 * pi

  ut <- rep(
    NA_real_,
    length(pt)
  )

  ut[
    is.finite(theta_all)
  ] <- predict(

    gam_mod,

    newdata = data.frame(
      theta =
        theta_all[
          is.finite(theta_all)
        ]
    )
  )

  rescale01(ut)
}

cat(
  "\nInitialisation terminée.\n"
)

############################################################
## DBSCAN
############################################################

run_dbscan <- function(
    emb,
    eps,
    minPts,
    keep_n = 2
){

  db <- dbscan::dbscan(
    emb,
    eps = eps,
    minPts = minPts
  )

  cl_tab <- table(db$cluster)

  cl_tab <- cl_tab[
    names(cl_tab) != "0"
  ]

  if(length(cl_tab) == 0){

    keep <- character(0)

  } else {

    keep <- names(
      sort(
        cl_tab,
        decreasing = TRUE
      )
    )[seq_len(
      min(
        keep_n,
        length(cl_tab)
      )
    )]

  }

  outlier <- !(
    as.character(db$cluster) %in% keep
  )

  list(
    cluster = db$cluster,
    outlier = outlier,
    keep = keep
  )
}

############################################################
## VISUALISATION OF OUTLIERS
############################################################

plot_dbscan_outliers <- function(
    obj,
    cluster_label,
    red_2d = "multimodal_umap_2D"
){

  cells <- colnames(obj)[
    obj$traj_cluster == cluster_label
  ]

  stopifnot(length(cells) > 0)

  flag_2d <- paste0(
    "out_dbscan__",
    red_2d,
    "__",
    cluster_label
  )

  stopifnot(
    flag_2d %in%
      colnames(obj@meta.data)
  )

  emb2 <- Embeddings(
    obj,
    reduction = red_2d
  )[cells, 1:2, drop = FALSE]

  df <- data.frame(
    x = emb2[, 1],
    y = emb2[, 2],
    outlier =
      obj@meta.data[
        cells,
        flag_2d
      ]
  )

  ggplot(
    df,
    aes(x, y)
  ) +

    geom_point(
      color = "grey85",
      size = 0.8,
      alpha = 0.8
    ) +

    geom_point(
      data = df[
        df$outlier,
        ,
        drop = FALSE
      ],
      color = "red",
      size = 1.2
    ) +

    theme_minimal(
      base_size = 12
    ) +

    labs(
      title = paste0(
        cluster_label,
        " — DBSCAN outliers"
      ),
      x = colnames(emb2)[1],
      y = colnames(emb2)[2]
    )
}


cat("\n====================\n")
cat("DBSCAN\n")
cat("====================\n")

for(cluster_label in clusters){

  cat("\n")
  cat(cluster_label, "\n")

  cells <- colnames(integrated73)[
    integrated73$traj_cluster ==
      cluster_label
  ]

  if(length(cells) < 50){

    cat("skip (<50 cellules)\n")
    next

  }

  emb2 <- Embeddings(
    integrated73,
    reduction = red_2d
  )[cells, 1:2, drop = FALSE]

  emb2 <- as.matrix(emb2)

  res2 <- run_dbscan(
    emb = emb2,
    eps = eps_2d,
    minPts = minPts,
    keep_n = keep_n_cores
  )

  cat(
    "keep cores :",
    paste(
      res2$keep,
      collapse = ", "
    ),
    "\n"
  )

  cat(
    "outliers :",
    sum(res2$outlier),
    "/",
    length(res2$outlier),
    "\n"
  )

  flag_2d <- paste0(
    "out_dbscan__",
    red_2d,
    "__",
    cluster_label
  )

  vec <- rep(
    FALSE,
    ncol(integrated73)
  )

  names(vec) <- colnames(
    integrated73
  )

  vec[cells] <- res2$outlier

  integrated73@meta.data[
    ,
    flag_2d
  ] <- vec

  ##########################################################
  ## EXPORT FIGURES
  ##########################################################

  p <- plot_dbscan_outliers(
    integrated73,
    cluster_label,
    red_2d
  )

  pdf(
    file.path(
      "results",
      "trajectory",
      "DBSCAN_outliers",
      paste0(
        cluster_label,
        "_DBSCAN_outliers.pdf"
      )
    ),
    width = 6,
    height = 5
  )

  print(p)

  dev.off()

  tiff(
    file.path(
      "results",
      "trajectory",
      "DBSCAN_outliers",
      paste0(
        cluster_label,
        "_DBSCAN_outliers.tiff"
      )
    ),
    width = 6,
    height = 5,
    units = "in",
    res = 300,
    compression = "lzw"
  )

  print(p)

  dev.off()
}

cat("\nDBSCAN terminé.\n")

############################################################
## VISUALISATION SLINGSHOT 2D
############################################################

plot_slingshot2d_pseudotime <- function(
    obj,
    cluster_label,
    red_2d = "multimodal_umap_2D"
){

  ##########################################################
  ## Cellules du cluster
  ##########################################################

  cells_all <- colnames(obj)[
    obj$traj_cluster == cluster_label
  ]

  stopifnot(length(cells_all) > 0)

  ##########################################################
  ## DBSCAN 2D
  ##########################################################

  flag_out <- paste0(
    "out_dbscan__",
    red_2d,
    "__",
    cluster_label
  )

  stopifnot(
    flag_out %in%
      colnames(obj@meta.data)
  )

  cells_keep <- cells_all[
    !obj@meta.data[
      cells_all,
      flag_out
    ]
  ]

  stopifnot(length(cells_keep) > 10)

  ##########################################################
  ## UMAP 2D
  ##########################################################

  umap2 <- Embeddings(
    obj,
    red_2d
  )[
    cells_keep,
    1:2,
    drop = FALSE
  ]

  ##########################################################
  ## SingleCellExperiment minimal
  ##########################################################

  sce2 <- SingleCellExperiment(
    assays = list(
      dummy = matrix(
        0,
        nrow = 1,
        ncol = length(cells_keep),
        dimnames = list(
          "dummy",
          cells_keep
        )
      )
    )
  )

  reducedDims(sce2)$UMAP2D <- umap2

  ##########################################################
  ## Slingshot
  ##########################################################

  sce2 <- slingshot(
    sce2,
    reducedDim = "UMAP2D",
    clusterLabels = rep(
      1,
      ncol(sce2)
    )
  )

  ##########################################################
  ## Pseudotime
  ##########################################################

  pt <- slingPseudotime(
    sce2
  )[,1]

  names(pt) <- colnames(
    sce2
  )

  ##########################################################
  ## Courbe Slingshot
  ##########################################################

  curve2 <- slingCurves(
    sce2
  )[[1]]

  curve2_df <- data.frame(
    x = curve2$s[
      curve2$ord,
      1
    ],
    y = curve2$s[
      curve2$ord,
      2
    ]
  )

  ##########################################################
  ## Cellules
  ##########################################################

  df <- data.frame(
    Cell = names(pt),
    x = umap2[
      names(pt),
      1
    ],
    y = umap2[
      names(pt),
      2
    ],
    pt = pt,
    stringsAsFactors = FALSE
  )

  df <- df[
    is.finite(df$pt),
    ,
    drop = FALSE
  ]

  ##########################################################
  ## Plot
  ##########################################################

  ggplot(
    df,
    aes(x, y)
  ) +

    geom_point(
      aes(color = pt),
      size = 1.2
    ) +

    scale_color_gradient(
      low = "blue",
      high = "yellow"
    ) +

    geom_path(
      data = curve2_df,
      aes(x, y),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 1.2
    ) +

    theme_minimal(
      base_size = 13
    ) +

    labs(
      title = paste0(
        "Slingshot 2D — ",
        cluster_label
      ),
      x = colnames(umap2)[1],
      y = colnames(umap2)[2],
      color = "pseudotime"
    )
}


############################################################
## ANCHORING ON D2
############################################################

anchor_pseudotime_on_D2 <- function(
    pt,
    cells,
    orig_ident,
    window_fraction = 0.05,
    d2_fraction = 0.10
){

  pt <- as.numeric(pt)
  names(pt) <- cells

 
  valid <- is.finite(pt)

  if(sum(valid) < 50){

    return(
      list(
        pt = pt,
        shift = NA_real_,
        window_size = NA_integer_,
        d2_threshold = NA_integer_,
        d2_in_window = NA_integer_,
        d2_proportion = NA_real_,
        anchor_position = NA_integer_,
        anchor_pseudotime = NA_real_,
        d2_start_cell = NA_character_
      )
    )
  }


  pt_valid <- pt[valid]
  cells_valid <- cells[valid]

  orig_valid <- orig_ident[cells_valid]

  ord <- order(pt_valid)

  ordered_cells <- cells_valid[ord]
  ordered_pt <- pt_valid[ord]
  ordered_orig <- orig_valid[ord]

  n <- length(ordered_cells)



  window_size <- ceiling(
    window_fraction * n
  )


  n_d2 <- sum(
    ordered_orig == "D2"
  )

  if(n_d2 == 0){

    cat(
      "Aucune cellule D2 : pas d'ancrage\n"
    )

    return(
      list(
        pt = pt,
        shift = 0,
        window_size = window_size,
        d2_threshold = NA_integer_,
        d2_in_window = 0,
        d2_proportion = NA_real_,
        anchor_position = NA_integer_,
        anchor_pseudotime = NA_real_,
        d2_start_cell = NA_character_
      )
    )
  }


  d2_threshold <- ceiling(
    d2_fraction * n_d2
  )


  is_d2 <- ordered_orig == "D2"

  anchor_position <- NA_integer_

  d2_in_window <- NA_integer_

  for(start in seq_len(n - window_size + 1)){

    end <- start + window_size - 1

    n_d2_window <- sum(
      is_d2[start:end]
    )

    if(n_d2_window >= d2_threshold){

      anchor_position <- start
      d2_in_window <- n_d2_window

      break
    }
  }


  if(is.na(anchor_position)){

    cat(
      "Aucune fenêtre contenant au moins",
      d2_threshold,
      "D2\n"
    )

    return(
      list(
        pt = pt,
        shift = 0,
        window_size = window_size,
        d2_threshold = d2_threshold,
        d2_in_window = NA_integer_,
        d2_proportion = NA_real_,
        anchor_position = NA_integer_,
        anchor_pseudotime = NA_real_,
        d2_start_cell = NA_character_
      )
    )
  }


  anchor_pseudotime <-
    ordered_pt[anchor_position]

  d2_start_cell <-
    ordered_cells[anchor_position]

  d2_proportion <-
    d2_in_window / window_size



  pt_norm <- rescale01(
    pt
  )

  ordered_pt_norm <- pt_norm[
    match(
      ordered_cells,
      cells
    )
  ]

  shift <- ordered_pt_norm[
    anchor_position
  ]

#

  pt_new <- pt_norm

  pt_new[valid] <- (
    pt_norm[valid] - shift
  ) %% 1

 

  pt_new[valid] <-
    pt_new[valid] * 100



  cat(
    "D2 total :",
    n_d2,
    "\n"
  )

  cat(
    "Seuil D2 :",
    d2_threshold,
    "\n"
  )

  cat(
    "Taille fenêtre :",
    window_size,
    "\n"
  )

  cat(
    "Première fenêtre valide :",
    anchor_position,
    "->",
    anchor_position + window_size - 1,
    "\n"
  )

  cat(
    "D2 dans fenêtre :",
    d2_in_window,
    "\n"
  )

  cat(
    "Proportion D2 :",
    round(
      100 * d2_proportion,
      2
    ),
    "%\n"
  )

  cat(
    "Pseudotime début :",
    anchor_pseudotime,
    "\n"
  )

  cat(
    "Cellule ancre :",
    d2_start_cell,
    "\n"
  )

  cat(
    "Shift :",
    round(
      100 * shift,
      3
    ),
    "\n"
  )



  list(
    pt = pt_new,
    shift = shift * 100,
    window_size = window_size,
    d2_threshold = d2_threshold,
    d2_in_window = d2_in_window,
    d2_proportion = d2_proportion,
    anchor_position = anchor_position,
    anchor_pseudotime = anchor_pseudotime,
    d2_start_cell = d2_start_cell
  )
}


############################################################
## SLINGSHOT 2D : pseudotime
############################################################

cat("\n====================\n")
cat("Slingshot 2D pseudotime\n")
cat("====================\n")

for(cluster_label in clusters){

  cat("\n")
  cat(cluster_label, "\n")

  flag_out <- paste0(
    "out_dbscan__",
    red_2d,
    "__",
    cluster_label
  )

  if(
    !(flag_out %in%
      colnames(
        integrated73@meta.data
      ))
  ){

    cat("DBSCAN absent\n")
    next

  }



  cells_all <- colnames(
    integrated73
  )[
    integrated73$traj_cluster ==
      cluster_label
  ]

  cells_keep <- cells_all[
    !integrated73@meta.data[
      cells_all,
      flag_out
    ]
  ]

  if(length(cells_keep) < 50){

    cat(
      "skip (<50 cellules après filtrage)\n"
    )

    next
  }



  umap2 <- Embeddings(
    integrated73,
    red_2d
  )[cells_keep, 1:2, drop = FALSE]


  sce2 <- SingleCellExperiment(
    assays = list(
      dummy = matrix(
        0,
        nrow = 1,
        ncol = length(cells_keep),
        dimnames = list(
          "dummy",
          cells_keep
        )
      )
    )
  )

  reducedDims(sce2)$UMAP2D <- umap2


  sce2 <- slingshot(
    sce2,
    reducedDim = "UMAP2D",
    clusterLabels = rep(
      1,
      ncol(sce2)
    )
  )

##########################################################
## Extraction pseudotime
##########################################################

pt_raw <- slingPseudotime(
  sce2
)[,1]

names(pt_raw) <- colnames(
  sce2
)



if(cluster_label == "Cluster_C"){

  pt_raw <- max(
    pt_raw,
    na.rm = TRUE
  ) - pt_raw

  cat(
    "Sens du pseudotemps inversé pour",
    cluster_label,
    "\n"
  )
}


col_pt_raw <- paste0(
  "pt_slingshot2d_raw_",
  cluster_label
)

integrated73@meta.data[
  cells_all,
  col_pt_raw
] <- NA_real_

integrated73@meta.data[
  names(pt_raw),
  col_pt_raw
] <- pt_raw


anchor <- anchor_pseudotime_on_D2(
  pt = pt_raw,
  cells = colnames(sce2),
  orig_ident = integrated73$orig.ident[
    colnames(sce2)
  ],
  window_fraction = 0.05,
  d2_fraction = 0.10
)

pt <- anchor$pt


cat(
  "D2 total :",
  sum(
    integrated73$orig.ident[
      names(pt_raw)
    ] == "D2"
  ),
  "\n"
)

cat(
  "Taille fenêtre :",
  anchor$window_size,
  "\n"
)

cat(
  "D2 dans fenêtre :",
  anchor$d2_in_window,
  "\n"
)

cat(
  "Proportion D2 :",
  round(
    100 * anchor$d2_proportion,
    2
  ),
  "%\n"
)

cat(
  "Position ancre :",
  anchor$anchor_position,
  "\n"
)

cat(
  "Pseudotime ancre :",
  anchor$anchor_pseudotime,
  "\n"
)

cat(
  "Cellule ancre :",
  anchor$d2_start_cell,
  "\n"
)

cat(
  "Shift circulaire :",
  anchor$shift,
  "\n"
)



pt <- pt[
  is.finite(pt)
]



ord <- order(
  pt,
  method = "radix"
)

ordered_cells <- names(pt)[ord]

n_ordered <- length(
  ordered_cells
)


pt_continuous <- seq(
  from = 0,
  to = 100,
  length.out = n_ordered
)

names(pt_continuous) <- ordered_cells



col_pt <- paste0(
  "pt_slingshot2d_",
  cluster_label
)

col_ord <- paste0(
  "ord_slingshot2d_",
  cluster_label
)

integrated73@meta.data[
  cells_all,
  col_pt
] <- NA_real_

integrated73@meta.data[
  cells_all,
  col_ord
] <- NA_integer_



integrated73@meta.data[
  ordered_cells,
  col_pt
] <- pt_continuous



integrated73@meta.data[
  ordered_cells,
  col_ord
] <- seq_len(
  n_ordered
)

cat(
  "Trajectoire continue sauvegardée :",
  n_ordered,
  "cellules\n"
)

cat(
  "Intervalle :",
  min(pt_continuous),
  "→",
  max(pt_continuous),
  "\n"
)


  cat(
    "pseudotime Slingshot 2D calculé pour",
    length(pt),
    "cellules\n"
  )



p <- tryCatch(

  plot_slingshot2d_pseudotime(
    integrated73,
    cluster_label,
    red_2d
  ),

  error = function(e){

    cat(
      "Erreur visualisation Slingshot 2D :",
      e$message,
      "\n"
    )

    NULL
  }
)

if(!is.null(p)){

  pdf(
    file.path(
      "results",
      "trajectory",
      "Slingshot2D",
      paste0(
        cluster_label,
        "_Slingshot2D_TrajectoryPlot.pdf"
      )
    ),
    width = 6,
    height = 5
  )

  print(p)

  dev.off()

  tiff(
    file.path(
      "results",
      "trajectory",
      "Slingshot2D",
      paste0(
        cluster_label,
        "_Slingshot2D_TrajectoryPlot.tiff"
      )
    ),
    width = 6,
    height = 5,
    units = "in",
    res = 300,
    compression = "lzw"
  )

  print(p)

  dev.off()
}

}

cat("\nPseudotime Slingshot 2D terminé.\n")



frise_slingshot2d <- list()

for(cl in clusters){

  cat(
    "Frise Slingshot 2D :",
    cl,
    "\n"
  )

  cells <- colnames(integrated73)[
    integrated73$traj_cluster == cl
  ]

  col_pt <- paste0(
    "pt_slingshot2d_",
    cl
  )

  if(
    !(col_pt %in%
      colnames(integrated73@meta.data))
  ){

    cat(
      "pseudotime absent\n"
    )

    next
  }

  df_frise <- tibble(
    Cell = cells,

    slingshot2d =
      integrated73@meta.data[
        cells,
        col_pt
      ],

    orig =
      integrated73$orig.ident[
        cells
      ]
  ) |>
    filter(
      is.finite(slingshot2d)
    ) |>
    arrange(
      slingshot2d
    )

  if(nrow(df_frise) < 30){

    cat(
      "skip (<30 cellules)\n"
    )

    next
  }

  frise_slingshot2d[[cl]] <-
    df_frise |>
    mutate(
      xmin =
        (row_number() - 1) / n(),

      xmax =
        row_number() / n(),

      cluster = cl
    )
}

frise_slingshot2d <- bind_rows(
  frise_slingshot2d
)

cat(
  "Cellules dans frise Slingshot 2D :",
  nrow(frise_slingshot2d),
  "\n"
)



############################################################
## ETA² SLINGSHOT 2D — PSEUDOTEMPS BRUT
############################################################

eta_slingshot2d <- list()

for(cl in clusters){

  cat(
    "ETA² Slingshot 2D brut :",
    cl,
    "\n"
  )

  cells <- colnames(integrated73)[
    integrated73$traj_cluster == cl
  ]


  col_pt_raw <- paste0(
    "pt_slingshot2d_raw_",
    cl
  )

  if(
    !(col_pt_raw %in%
      colnames(integrated73@meta.data))
  ){

    cat(
      "pseudotime Slingshot brut absent\n"
    )

    next
  }

  df_eta <- tibble(
    Cell = cells,

    slingshot2d_raw =
      integrated73@meta.data[
        cells,
        col_pt_raw
      ],

    orig =
      integrated73$orig.ident[
        cells
      ]
  ) |>
    filter(
      is.finite(slingshot2d_raw)
    )

  if(nrow(df_eta) < 30){

    cat(
      "skip (<30 cellules)\n"
    )

    next
  }

  eta_slingshot2d[[cl]] <-
    eta2_from_aov(
      df_eta$slingshot2d_raw,
      df_eta$orig
    )

  cat(
    "ETA² brut =",
    eta_slingshot2d[[cl]],
    "\n"
  )
}


cat("\nSauvegarde objet final...\n")

saveRDS(
  integrated73,
  file =
    "results/objects/PGintegrated73.reclustered.linkpeaks.trajectory.rds"
)

cat(
  "Objet sauvegardé :\n",
  "results/objects/PGintegrated73.reclustered.linkpeaks.trajectory.rds\n"
)

cat("\n07_trajectory.R terminé.\n")