
suppressPackageStartupMessages({

  library(data.table)
  library(pheatmap)
  library(viridis)
  library(segmented)

})



q_cutoff <- 0.05

k_force_RNA <- c(
  Cluster_G   = 4
)

k_force_ATAC <- c(
  Cluster_G   = NA
)

k_force_Velocity <- c(
  Cluster_G   = NA
)

spline_spar <- 0.65

k_max <- 8



dir.create(
  "results/temporal_modules",
  recursive = TRUE,
  showWarnings = FALSE
)



minmax_11 <- function(x){

  rng <- range(x)

  if(diff(rng) == 0)
    return(rep(0,length(x)))

  2 * (x - rng[1]) / diff(rng) - 1

}

############################################################
## PREPARE MATRIX
############################################################

prepare_heatmap_matrix <- function(
    mat,
    spline_spar = 0.45
){

  time_bins <- seq_len(
    ncol(mat)
  )

  mat_smooth <- t(

    apply(

      mat,

      1,

      function(e){

        smooth.spline(

          x = time_bins,
          y = e,
          spar = spline_spar

        )$y

      }

    )

  )

  rownames(mat_smooth) <- rownames(mat)

  mat_scaled <- t(
    scale(
      t(mat_smooth)
    )
  )

  mat_scaled[
    is.na(mat_scaled)
  ] <- 0

  mat_scaled <- t(

    apply(

      mat_scaled,

      1,

      minmax_11

    )

  )

  rownames(mat_scaled) <- rownames(mat)

  colnames(mat_scaled) <- colnames(mat)

  mat_scaled

}



segm_heatmap <- function(
    mat,
    obj_label,
    k_force = NA,
    heat_cols,
    k_max = 8
){

  min_len <- min(
  20,
  floor(nrow(mat) / 3)
)

min_len <- max(
  min_len,
  5
)

  ##########################################################
  ## PLATEAU
  ##########################################################

  plateau <- function(
      x,
      thr = 0.85,
      min_run = 2
  ){

    m <- max(x)

    if(m == 0)
      return(which.max(x))

    flag <- x/m >= thr

    r <- rle(flag)

    idx <- which(
      r$values
    )

    if(!length(idx))
      return(which.max(x))

    st <- cumsum(
      c(
        1,
        head(
          r$lengths,
          -1
        )
      )
    )[idx]

    en <- st +
      r$lengths[idx] - 1

    ok <- which(
      (en-st+1) >= min_run
    )

    if(!length(ok))
      return(which.max(x))

    mid <- (
      st[ok] +
      en[ok]
    ) / 2

    mid[
      which.max(
        r$lengths[idx][ok]
      )
    ]

  }

  ##########################################################
  ## ORDRE
  ##########################################################

  moment <- apply(
    mat,
    1,
    plateau
  )

  ord <- names(
    sort(moment)
  )

  y <- moment[ord]

  x <- seq_along(ord)

  lm0 <- lm(
    y ~ x
  )

 #########################################################
  ## BIC
########################################################

cat("Calcul BIC\n")
flush.console()

models <- vector(
  "list",
  k_max + 1
)

bic <- sapply(
  0:k_max,
  function(k){

    mdl <- if(k == 0){

      lm0

    } else {

      tryCatch(

        segmented(
          lm0,
          seg.Z = ~x,
          npsi = k,
          control = seg.control(
            display = FALSE,
            it.max = 50
          )
        ),

        error = function(e) lm0,
        warning = function(w) lm0

      )

    }

    models[[k + 1]] <<- mdl

    BIC(mdl)

  }
)

print(
  data.frame(
    k = 0:k_max,
    BIC = round(bic, 2),
    delta_best = round(
      bic - min(bic),
      2
    ),
    delta_step = round(
      c(NA, diff(bic)),
      2
    )
  )
)

flush.console()

cat("BIC OK\n")
flush.console()



 best_k_auto <- which.min(bic)-1

cat(
  "best_k_auto =",
  best_k_auto,
  "\n"
)

if(!is.na(k_force))
  best_k <- k_force
else
  best_k <- best_k_auto

if(is.na(k_force))
  best_k <- min(best_k, 5)

cat(
  "best_k_final =",
  best_k,
  "\n"
)

flush.console()

  ##########################################################
  ## BREAKS
  ##########################################################

  breaks <- integer(0)

cat("ENTER BREAKS\n")
flush.console()

  if(best_k > 0){

  mdl_best <- models[[best_k + 1]]

  if(
    inherits(mdl_best, "segmented") &&
    !is.null(mdl_best$psi)
  ){

    breaks <- sort(
      mdl_best$psi[, "Est."]
    )

  }

}



  ##########################################################
  ## SEGMENTS
  ##########################################################

  seg_id <- cut(

    x,

    c(
      0,
      breaks,
      length(ord)
    ),

    labels = FALSE,

    include.lowest = TRUE

  )

  K <- length(
    unique(seg_id)
  )

  ann <- data.frame(
    Module = factor(seg_id)
  )

  rownames(ann) <- ord

  ##########################################################
  ## HEATMAP
  ##########################################################

module_levels <- sort(
  unique(seg_id)
)

module_cols <- viridis::viridis(
  length(module_levels)
)

names(module_cols) <- module_levels

ann_colors <- list(
  Module = module_cols
)

  hm <- pheatmap(

    mat[ord,,drop=FALSE],

    cluster_rows = FALSE,
    cluster_cols = FALSE,

annotation_row = ann,

annotation_colors = ann_colors,

    show_rownames = FALSE,
    show_colnames = FALSE,

    color = heat_cols,

    main = paste0(
      obj_label,
      " (",
      K,
      " modules)"
    )

  )

  list(

    heatmap = hm,

    ordered_matrix =
      mat[ord,,drop=FALSE],

    seg_id = seg_id,

    genes_order = ord,

    K = K,

    BIC = bic,

    breaks = breaks

  )

}



DE_results <- readRDS(
  "results/pseudotime_DE/objects/DE_results.rds"
)

Velocity_results <- readRDS(
  "results/pseudotime_DE/objects/Velocity_results.rds"
)

clusters <- names(
  DE_results
)





for(cl in clusters){

  cat(
    "\n====================\n"
  )

  cat(
    "Cluster:",
    cl,
    "\n"
  )

  cat(
    "====================\n"
  )

  outdir <- file.path(
    "results/temporal_modules",
    cl
  )

  dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  ##########################################################
  ## VALIDATED RNA
  ##########################################################

  RNA_perm <- fread(

    file.path(

      "results/pseudotime_perm/RNA",

      paste0(
        cl,
        "_RNA_perm.csv"
      )

    )

  )

  RNA_valid <- RNA_perm[
    q_emp < q_cutoff
  ]$gene


expr_keep <-
    DE_results[[cl]]$expr_keep

  expr_keep <-
    expr_keep[
      intersect(
        rownames(expr_keep),
        RNA_valid
      ),
      ,
      drop = FALSE
    ]

cat("RNA valid:", nrow(expr_keep), "\n")
flush.console()

##########################################################
## VALIDATED VELOCITY
##########################################################

Velocity_file <- file.path(
  "results/pseudotime_perm/Velocity",
  paste0(
    cl,
    "_Velocity_perm.csv"
  )
)

if(file.exists(Velocity_file)){

  Velocity_perm <- fread(
    Velocity_file
  )

  Velocity_valid <- Velocity_perm[
    q_emp < q_cutoff
  ]$gene

  velocity_keep <- Velocity_results[[cl]]$velocity_filt

  velocity_keep <- velocity_keep[
    intersect(
      rownames(velocity_keep),
      Velocity_valid
    ),
    ,
    drop = FALSE
  ]

  cat(
    "Velocity valid:",
    nrow(velocity_keep),
    "\n"
  )

} else {

  velocity_keep <- NULL

  cat(
    "Velocity permutation absente — Velocity ignorée\n"
  )

}

flush.console()

##########################################################
## RNA HEATMAP
##########################################################

if(nrow(expr_keep) == 0){

  cat(
    "No validated RNA genes for",
    cl,
    "\n"
  )

  flush.console()

} else {

  cat("Prepare RNA\n")
  flush.console()

  RNA_heatmap <- prepare_heatmap_matrix(
    expr_keep,
    spline_spar
  )

  cat("RNA ready\n")
  flush.console()

  cat("Segment RNA\n")
  flush.console()

  k_rna <- k_force_RNA[[cl]]

  cat(
    "RNA k_force =",
    k_rna,
    "\n"
  )

  flush.console()

  RNA_res <- segm_heatmap(

    RNA_heatmap,

    paste0(
      cl,
      " RNA"
    ),

    k_force = k_rna,

    heat_cols =
      viridis::plasma(50),

    k_max = k_max

  )

  cat("RNA segmented\n")
  flush.console()

  tiff(
    file.path(
      outdir,
      "RNA_heatmap.tiff"
    ),
    width = 10,
    height = 8,
    units = "in",
    res = 300,
    compression = "lzw"
  )

  print(
    RNA_res$heatmap
  )

  dev.off()

  RNA_modules <- data.frame(

    gene =
      RNA_res$genes_order,

    module =
      RNA_res$seg_id

  )

  fwrite(
    RNA_modules,
    file.path(
      outdir,
      "RNA_gene_modules.csv"
    )
  )

  saveRDS(
    RNA_res,
    file.path(
      outdir,
      "RNA_segmented.rds"
    )
  )
}


##########################################################
## VALIDATED ATAC
##########################################################

Peak_file <- file.path(
  "results/pseudotime_perm/Peaks",
  paste0(
    cl,
    "_Peak_perm.csv"
  )
)

if(file.exists(Peak_file)){

  Peak_perm <- fread(
    Peak_file
  )

  Peak_valid <- Peak_perm[
    q_emp < q_cutoff
  ]$peak

  peak_keep <- DE_results[[cl]]$peak_keep

  peak_keep <- peak_keep[
    intersect(
      rownames(peak_keep),
      Peak_valid
    ),
    ,
    drop = FALSE
  ]

  cat(
    "ATAC valid:",
    nrow(peak_keep),
    "\n"
  )

} else {

  peak_keep <- NULL

  cat(
    "ATAC permutation absente — ATAC ignoré\n"
  )

}

flush.console()


############################################################
## ATAC HEATMAP — VALIDATED PEAKS
############################################################

if(
  is.null(peak_keep) ||
  nrow(peak_keep) == 0
){

  cat(
    "No validated ATAC peaks for",
    cl,
    "\n"
  )

  flush.console()

} else {

  cat(
    "Prepare ATAC heatmap\n"
  )

  cat(
    "Peaks utilisés :",
    nrow(peak_keep),
    "\n"
  )

  flush.console()


  ##########################################################
  ## SMOOTHING + NORMALISATION
  ##########################################################

  Peak_heatmap <- prepare_heatmap_matrix(
    peak_keep,
    spline_spar
  )

  cat(
    "ATAC heatmap matrix :",
    nrow(Peak_heatmap),
    "peaks x",
    ncol(Peak_heatmap),
    "bins\n"
  )

  flush.console()


  ##########################################################
  ## SEGMENTATION
  ##########################################################

  k_atac <- k_force_ATAC[[cl]]

  cat(
    "Segment ATAC\n",
    "k_force =",
    k_atac,
    "\n"
  )

  flush.console()

  Peak_res <- segm_heatmap(

    Peak_heatmap,

    paste0(
      cl,
      " ATAC"
    ),

    k_force = k_atac,

    heat_cols =
      viridis::mako(50),

    k_max = k_max

  )

  cat(
    "ATAC segmented\n",
    "Modules =",
    Peak_res$K,
    "\n"
  )

  flush.console()


  ##########################################################
  ## TABLE PEAK → MODULE
  ##########################################################

  Peak_modules <- data.table(

    peak =
      Peak_res$genes_order,

    module =
      Peak_res$seg_id,

    order =
      seq_along(
        Peak_res$genes_order
      )

  )


  ##########################################################
  ## SAUVEGARDE MODULES
  ##########################################################

  fwrite(

    Peak_modules,

    file.path(
      outdir,
      "Peak_modules.csv"
    )

  )


  ##########################################################
  ## SAUVEGARDE MATRICE ORDONNÉE
  ##########################################################

  saveRDS(

    Peak_res$ordered_matrix,

    file.path(
      outdir,
      "Peak_heatmap_matrix.rds"
    )

  )


  ##########################################################
  ## SAUVEGARDE OBJET COMPLET
  ##########################################################

  saveRDS(

    Peak_res,

    file.path(
      outdir,
      "Peak_segmented.rds"
    )

  )


  ##########################################################
  ## TIFF
  ##########################################################

  tiff(

    file.path(
      outdir,
      "Peak_heatmap.tiff"
    ),

    width = 12,
    height = 10,

    units = "in",

    res = 300,

    compression = "lzw"

  )

  grid::grid.newpage()

  grid::grid.draw(
    Peak_res$heatmap$gtable
  )

  dev.off()


  ##########################################################
  ## BIC
  ##########################################################

  fwrite(

    data.table(

      k =
        0:length(Peak_res$BIC) - 1,

      BIC =
        Peak_res$BIC

    ),

    file.path(
      outdir,
      "Peak_BIC.csv"
    )

  )


  cat(
    "ATAC heatmap saved\n"
  )

  flush.console()

}

##########################################################
## VELOCITY / GDV HEATMAP
##########################################################

if(
  is.null(velocity_keep) ||
  nrow(velocity_keep) == 0
){
  cat(
    "No validated Velocity genes for",
    cl,
    "\n"
  )

  flush.console()

} else {

  cat("Prepare Velocity\n")
  flush.console()

  Velocity_heatmap <- prepare_heatmap_matrix(
    velocity_keep,
    spline_spar
  )

  cat("Velocity ready\n")
  flush.console()

  cat("Segment Velocity\n")
  flush.console()

k_velocity <- k_force_Velocity[[cl]]

  cat(
    "Velocity k_force =",
    k_velocity,
    "\n"
  )

  flush.console()

  Velocity_res <- segm_heatmap(

    Velocity_heatmap,

    paste0(
      cl,
      " Velocity"
    ),

    k_force = k_velocity,

    heat_cols =
      viridis::viridis(50),

    k_max = k_max

  )

  cat("Velocity segmented\n")
  flush.console()

  tiff(
    file.path(
      outdir,
      "Velocity_heatmap.tiff"
    ),
    width = 10,
    height = 8,
    units = "in",
    res = 300,
    compression = "lzw"
  )

  print(
    Velocity_res$heatmap
  )

  dev.off()

  Velocity_modules <- data.frame(

    gene =
      Velocity_res$genes_order,

    module =
      Velocity_res$seg_id

  )

  fwrite(
    Velocity_modules,
    file.path(
      outdir,
      "Velocity_gene_modules.csv"
    )
  )

  saveRDS(
    Velocity_res,
    file.path(
      outdir,
      "Velocity_segmented.rds"
    )
  )
}

}