
suppressPackageStartupMessages({

  library(Seurat)
  library(Signac)
  library(Matrix)
  library(data.table)

})

options(
  stringsAsFactors = FALSE
)


trajectory_prefix <-
  "pt_slingshot2d_"

cells_per_bin <-
  10

gam_k <-
  5

gam_qvalue <-
  0.05



nperm_task <-
  125

nperm_total <-
  1000


task_id <-
  as.integer(
    Sys.getenv(
      "SLURM_ARRAY_TASK_ID",
      "1"
    )
  )

n_tasks <-
  as.integer(
    Sys.getenv(
      "SLURM_ARRAY_TASK_COUNT",
      "8"
    )
  )

if(
  task_id < 1 ||
  task_id > n_tasks
){

  stop(
    "SLURM_ARRAY_TASK_ID invalide : ",
    task_id
  )

}


seed <-
  20260814 +
  task_id

set.seed(
  seed
)


out_dir <-
  "results/pseudotime_perm/RNA"

object_dir <-
  "results/pseudotime_perm/objects/RNA"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  object_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


cat(
  "\n========================================\n",
  "09 — PERMUTATIONS RNA\n",
  "========================================\n"
)

cat(
  "Tâche SLURM : ",
  task_id,
  " / ",
  n_tasks,
  "\n",
  sep = ""
)

cat(
  "Permutations cette tâche : ",
  nperm_task,
  "\n",
  sep = ""
)

cat(
  "Permutations totales attendues : ",
  nperm_total,
  "\n",
  sep = ""
)

cat(
  "Seed : ",
  seed,
  "\n",
  sep = ""
)


cat(
  "\n====================\n",
  "Chargement\n",
  "====================\n"
)

obj <-
  readRDS(
    "results/objects/PGintegrated73.reclustered.linkpeaks.trajectory.rds"
  )

DE_results <-
  readRDS(
    "results/pseudotime_DE/objects/DE_results.rds"
  )

cat(
  "Objet chargé\n"
)

clusters <-
  c(
    "Cluster_L"
  )

if(
  length(clusters) != 1
){

  stop(
    "Ce script doit traiter un seul cluster."
  )

}

cl <-
  clusters[1]

cat(
  "Cluster : ",
  cl,
  "\n",
  sep = ""
)


if(
  !cl %in% names(DE_results)
){

  stop(
    "Cluster absent de DE_results : ",
    cl
  )

}


bin_expression <-
  function(
    mat,
    traj,
    cells_per_bin = 10
  ){

    ord <-
      order(
        traj
      )

    mat <-
      mat[
        ,
        ord,
        drop = FALSE
      ]

    traj <-
      traj[
        ord
      ]

    bins <-
      split(
        seq_along(traj),
        ceiling(
          seq_along(traj) /
          cells_per_bin
        )
      )

    binned <-
      sapply(
        bins,
        function(idx){

          Matrix::rowMeans(
            mat[
              ,
              idx,
              drop = FALSE
            ]
          )

        }
      )

    binned <-
      as.matrix(
        binned
      )

    ## Cas particulier : un seul bin

    if(
      is.null(
        dim(binned)
      )
    ){

      binned <-
        matrix(
          binned,
          nrow = nrow(mat)
        )

    }

    rownames(binned) <-
      rownames(mat)

    binned

  }


compute_GAMstats <-
  function(
    mat,
    traj
  ){


    expr_bin <-
      bin_expression(
        mat =
          mat,
        traj =
          traj,
        cells_per_bin =
          cells_per_bin
      )

    time_bins <-
      seq_len(
        ncol(expr_bin)
      )


    result_R2 <-
      numeric(
        nrow(expr_bin)
      )

    result_dev <-
      numeric(
        nrow(expr_bin)
      )

    names(result_R2) <-
      rownames(
        expr_bin
      )

    names(result_dev) <-
      rownames(
        expr_bin
      )


    for(
      i in seq_len(
        nrow(expr_bin)
      )
    ){

      e <-
        as.numeric(
          expr_bin[
            i,
            ]
        )

      if(
        anyNA(e) ||
        !is.finite(
          sd(e)
        ) ||
        sd(e) == 0
      ){

        result_R2[i] <-
          0

        result_dev[i] <-
          0

        next

      }

      df <-
        data.frame(

          expression =
            e,

          pseudotime =
            time_bins

        )

      fit <-
        try(

          mgcv::gam(

            expression ~
              s(
                pseudotime,
                k = gam_k,
                bs = "tp"
              ),

            data =
              df,

            method =
              "REML"

          ),

          silent = TRUE

        )

      if(
        inherits(
          fit,
          "try-error"
        )
      ){

        result_R2[i] <-
          0

        result_dev[i] <-
          0

        next

      }

      sm <-
        summary(
          fit
        )

      r2 <-
        sm$r.sq

      dev <-
        sm$dev.expl

      if(
        !is.finite(r2)
      ){

        r2 <-
          0

      }

      if(
        !is.finite(dev)
      ){

        dev <-
          0

      }

      result_R2[i] <-
        r2

      result_dev[i] <-
        dev

    }

    list(

      R2 =
        result_R2,

      deviance =
        result_dev

    )

  }



cat(
  "\n====================\n",
  "Préparation données\n",
  "====================\n"
)


cells <-
  colnames(obj)[
    obj$traj_cluster == cl
  ]

cat(
  "Cellules cluster : ",
  length(cells),
  "\n",
  sep = ""
)


traj_col <-
  paste0(
    trajectory_prefix,
    cl
  )

if(
  !traj_col %in%
  colnames(obj@meta.data)
){

  stop(
    "Trajectoire absente : ",
    traj_col
  )

}

traj <-
  obj@meta.data[
    cells,
    traj_col
  ]

names(traj) <-
  cells


keep <-
  is.finite(
    traj
  )

cells <-
  cells[
    keep
  ]

traj <-
  traj[
    keep
  ]

cat(
  "Cellules avec pseudotemps : ",
  length(traj),
  "\n",
  sep = ""
)


traj <-
  traj -
  min(
    traj,
    na.rm = TRUE
  )

traj <-
  traj /
  max(
    traj,
    na.rm = TRUE
  ) *
  100

cat(
  "\n====================\n",
  "RNA\n",
  "====================\n"
)


RNA_genes <-
  DE_results[[cl]]$RNA_dyn

RNA_genes <-
  intersect(

    RNA_genes,

    rownames(
      obj[["SoupXRNA"]]
    )

  )

cat(
  "RNA dyn : ",
  length(RNA_genes),
  "\n",
  sep = ""
)

if(
  length(RNA_genes) == 0
){

  stop(
    "Aucun RNA_dyn pour ",
    cl
  )

}


DefaultAssay(
  obj
) <-
  "SoupXRNA"

expr_cell <-
  LayerData(

    object =
      obj,

    assay =
      "SoupXRNA",

    layer =
      "data"

  )[

    RNA_genes,
    cells,
    drop = FALSE

  ]

if(
  inherits(
    expr_cell,
    "dgCMatrix"
  )
){

  expr_cell <-
    as.matrix(
      expr_cell
    )

}

cat(
  "RNA : ",
  nrow(expr_cell),
  " genes x ",
  ncol(expr_cell),
  " cellules\n",
  sep = ""
)


cat(
  "\n====================\n",
  "GAM OBSERVE\n",
  "====================\n"
)

cat(
  "Calcul GAM observé RNA\n"
)

GAM_real <-
  compute_GAMstats(
    expr_cell,
    traj
  )

cat(
  "RNA OK\n"
)

GAM_real_R2 <-
  GAM_real$R2

GAM_real_dev <-
  GAM_real$deviance


RNA_observed <-
  data.table(

    gene =
      names(
        GAM_real_R2
      ),

    R2_obs =
      as.numeric(
        GAM_real_R2
      ),

    deviance_obs =
      as.numeric(
        GAM_real_dev
      )

  )


cat(
  "\n====================\n",
  "PERMUTATIONS\n",
  "====================\n"
)

cat(
  "Nombre total : ",
  nperm_total,
  "\n",
  sep = ""
)

cat(
  "Nombre cette tâche : ",
  nperm_task,
  "\n",
  sep = ""
)


perm_R2 <-
  matrix(

    NA_real_,

    nrow =
      length(
        GAM_real_R2
      ),

    ncol =
      nperm_task,

    dimnames =
      list(

        names(
          GAM_real_R2
        ),

        paste0(
          "perm_",
          seq_len(
            nperm_task
          )
        )

      )

  )

perm_deviance <-
  matrix(

    NA_real_,

    nrow =
      length(
        GAM_real_dev
      ),

    ncol =
      nperm_task,

    dimnames =
      list(

        names(
          GAM_real_dev
        ),

        paste0(
          "perm_",
          seq_len(
            nperm_task
          )
        )

      )

  )

for(
  i in seq_len(
    nperm_task
  )
){

  if(
    i == 1 ||
    i %% 25 == 0 ||
    i == nperm_task
  ){

    cat(
      "Permutation ",
      i,
      " / ",
      nperm_task,
      "\n",
      sep = ""
    )

  }

  traj_perm <-
    sample(
      traj,
      replace = FALSE
    )

  GAM_perm <-
    compute_GAMstats(
      expr_cell,
      traj_perm
    )



  perm_R2[
    ,
    i
  ] <-
    GAM_perm$R2



  perm_deviance[
    ,
    i
  ] <-
    GAM_perm$deviance

}



cat(
  "\nDimensions R² permutation :\n"
)

print(
  dim(
    perm_R2
  )
)

cat(
  "\nDimensions déviance permutation :\n"
)

print(
  dim(
    perm_deviance
  )
)


R2_perm_median <-
  apply(
    perm_R2,
    1,
    median,
    na.rm = TRUE
  )

R2_perm_q95 <-
  apply(
    perm_R2,
    1,
    quantile,
    probs = 0.95,
    na.rm = TRUE
  )

R2_perm_q99 <-
  apply(
    perm_R2,
    1,
    quantile,
    probs = 0.99,
    na.rm = TRUE
  )

dev_perm_median <-
  apply(
    perm_deviance,
    1,
    median,
    na.rm = TRUE
  )

dev_perm_q95 <-
  apply(
    perm_deviance,
    1,
    quantile,
    probs = 0.95,
    na.rm = TRUE
  )


RNA_chunk_summary <-
  data.table(

    gene =
      names(
        GAM_real_R2
      ),

    R2_obs =
      as.numeric(
        GAM_real_R2
      ),

    deviance_obs =
      as.numeric(
        GAM_real_dev
      ),

    R2_perm_median =
      as.numeric(
        R2_perm_median
      ),

    R2_perm_q95 =
      as.numeric(
        R2_perm_q95
      ),

    R2_perm_q99 =
      as.numeric(
        R2_perm_q99
      ),

    delta_R2 =
      as.numeric(
        GAM_real_R2 -
        R2_perm_median
      ),

    R2_enrichment_q95 =
      as.numeric(
        GAM_real_R2 /
        pmax(
          R2_perm_q95,
          .Machine$double.eps
        )
      ),

    dev_perm_median =
      as.numeric(
        dev_perm_median
      ),

    dev_perm_q95 =
      as.numeric(
        dev_perm_q95
      )

  )


############################################################
## EXPORT CSV CHUNK
############################################################

csv_file <-
  file.path(

    out_dir,

    paste0(
      cl,
      "_task",
      task_id,
      "_RNA_permutations.csv"
    )

  )

fwrite(
  RNA_chunk_summary,
  csv_file
)


rds_file <-
  file.path(

    object_dir,

    paste0(
      cl,
      "_task",
      task_id,
      "_RNA_permutations.rds"
    )

  )

saveRDS(

  list(

    cluster =
      cl,

    task_id =
      task_id,

    n_tasks =
      n_tasks,

    nperm_task =
      nperm_task,

    nperm_total =
      nperm_total,

    seed =
      seed,

    trajectory =
      traj_col,

    cells_per_bin =
      cells_per_bin,

    gam_k =
      gam_k,

    statistic =
      "R2_and_deviance_explained",

    observed =
      RNA_observed,

    perm_R2 =
      perm_R2,

    perm_deviance =
      perm_deviance

  ),

  rds_file

)



cat(
  "\n========================================\n",
  "TACHE TERMINEE\n",
  "========================================\n"
)

cat(
  "Cluster : ",
  cl,
  "\n",
  sep = ""
)

cat(
  "Task : ",
  task_id,
  "\n",
  sep = ""
)

cat(
  "Permutations : ",
  nperm_task,
  "\n",
  sep = ""
)

cat(
  "Fichier CSV : ",
  csv_file,
  "\n",
  sep = ""
)

cat(
  "Fichier RDS : ",
  rds_file,
  "\n",
  sep = ""
)

cat(
  "========================================\n"
)