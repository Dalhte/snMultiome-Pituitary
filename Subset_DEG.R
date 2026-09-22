suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(Matrix)
})

set.seed(1234)

############################################################
# PARAMETERS
############################################################

INPUT <- "/shared/projects/femcycle/PG_multiome_cluster/results/PGintegrated73.reclustered.rds"

OUTDIR <- "/shared/projects/femcycle/PG_multiome_cluster/results/deseq2_stage"

# Metadata
STAGE_COL  <- "orig.ident"
SAMPLE_COL <- "sample"
CLUSTER_COL <- "cluster_pb"

# Minimum number of cells required for a
# sample x cluster pseudobulk
MIN_CELLS <- 1

# Minimum number of biological replicates
# required in EACH condition
MIN_REPLICATES <- 2

# Statistical thresholds
FDR_THRESHOLD <- 0.05

LFC_THRESHOLD <- 1


# Clusters to subsample to the per-sample cell number of Cluster_G
SUBSAMPLE_CLUSTERS <- c(
  "Cluster_FSC",
  "Cluster_L",
  "Cluster_S"
)

REFERENCE_CLUSTER <- "Cluster_G"

# Only these comparisons are performed
COMPARISONS <- list(
  c("D2_PM", "PS"),
  c("PS", "E"),
  c("E", "D2_PM")
)

# Analyses performed sequentially
ASSAYS <- c(
  "RNA"
)

############################################################
# LOAD
############################################################

cat("========================================================\n")
cat("03.5 DESeq2 STAGE ANALYSIS\n")
cat("========================================================\n")

cat("Loading object...\n")

obj <- readRDS(INPUT)

cat(
  "Cells:",
  ncol(obj),
  "\n"
)

cat(
  "Features RNA:",
  nrow(obj[["SoupXRNA"]]),
  "\n"
)

cat(
  "Features ATAC:",
  nrow(obj[["peaks"]]),
  "\n"
)

############################################################
# METADATA
############################################################

meta <- obj@meta.data

required.metadata <- c(
  STAGE_COL,
  SAMPLE_COL,
  "new_cluster"
)

missing.metadata <- setdiff(
  required.metadata,
  colnames(meta)
)

if (length(missing.metadata) > 0) {

  stop(
    "Missing metadata columns: ",
    paste(
      missing.metadata,
      collapse = ", "
    )
  )
}

############################################################
# CLUSTER COLLAPSING
############################################################

meta[[CLUSTER_COL]] <- as.character(
  meta$new_cluster
)

meta[[CLUSTER_COL]][
  meta[[CLUSTER_COL]] %in%
    c("Cluster_G1", "Cluster_G2")
] <- "Cluster_G"

meta[[CLUSTER_COL]][
  meta[[CLUSTER_COL]] %in%
    c("Cluster_L1", "Cluster_L2")
] <- "Cluster_L"

obj$cluster_pb <- meta[[CLUSTER_COL]]

meta <- obj@meta.data

############################################################
# CHECK STAGES
############################################################

cat("\n========================================================\n")
cat("STAGES\n")
cat("========================================================\n")

cat(
  "Observed stages:\n"
)

print(
  table(
    meta[[STAGE_COL]],
    useNA = "ifany"
  )
)

requested.stages <- unique(
  unlist(COMPARISONS)
)

observed.stages <- unique(
  as.character(
    meta[[STAGE_COL]]
  )
)

observed.stages_collapsed <- unique(
  ifelse(
    observed.stages %in% c("D2", "PM"),
    "D2_PM",
    observed.stages
  )
)

missing.stages <- setdiff(
  requested.stages,
  observed.stages_collapsed
)

if (length(missing.stages) > 0) {

  stop(
    "Missing requested stages: ",
    paste(
      missing.stages,
      collapse = ", "
    )
  )
}

############################################################
# CHECK SAMPLE -> ORIG.IDENT
############################################################

sample.stage.table <- unique(
  meta[
    ,
    c(
      SAMPLE_COL,
      STAGE_COL
    ),
    drop = FALSE
  ]
)

sample.stage.counts <- table(
  sample.stage.table[[SAMPLE_COL]]
)

bad.samples <- names(
  sample.stage.counts[
    sample.stage.counts > 1
  ]
)

if (length(bad.samples) > 0) {

  stop(
    "ERROR: the same sample is associated with ",
    "more than one orig.ident:\n",
    paste(
      bad.samples,
      collapse = ", "
    )
  )
}

############################################################
# OUTPUT DIRECTORIES
############################################################

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(
    OUTDIR,
    "RNA"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(
    OUTDIR,
    "ATAC"
  ),
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
# FUNCTION
############################################################

run_stage_analysis <- function(
  object,
  assay_name
) {

  cat("\n\n")
  cat("########################################################\n")
  cat("ASSAY:", assay_name, "\n")
  cat("########################################################\n")

  DefaultAssay(object) <- assay_name

  counts <- GetAssayData(
    object,
    assay = assay_name,
    layer = "counts"
  )

  meta <- object@meta.data

  clusters <- sort(
    unique(
      meta[[CLUSTER_COL]]
    )
  )
  ##########################################################
  # REFERENCE CELL COUNTS FROM CLUSTER_G
  #
  # These counts define the subsampling target for
  # Cluster_FSC, Cluster_L and Cluster_S.
  ##########################################################

  reference_cells <- table(
    meta[
      meta[[CLUSTER_COL]] == REFERENCE_CLUSTER,
      SAMPLE_COL
    ]
  )

  reference_cells <- as.integer(
    reference_cells
  )

  names(reference_cells) <- names(
    table(
      meta[
        meta[[CLUSTER_COL]] == REFERENCE_CLUSTER,
        SAMPLE_COL
      ]
    )
  )

  cat("\n")
  cat("========================================================\n")
  cat(
    "REFERENCE CELL COUNTS:",
    REFERENCE_CLUSTER,
    "\n"
  )
  cat("========================================================\n")

  print(
    reference_cells
  )
  assay_outdir <- file.path(
    OUTDIR,
    ifelse(
      assay_name == "SoupXRNA",
      "RNA",
      "ATAC"
    )
  )

  ##########################################################
  # STORAGE
  ##########################################################

  comparison_summary <- list()

  unique_features_by_cluster <- list()

  ##########################################################
  # CLUSTERS
  ##########################################################

  for (cl in clusters) {

    cat("\n")
    cat("========================================================\n")
    cat("CLUSTER:", cl, "\n")
    cat("========================================================\n")

    cluster_cells <- rownames(meta)[
      meta[[CLUSTER_COL]] == cl
    ]

    cat(
      "Total cells:",
      length(cluster_cells),
      "\n"
    )

    if (length(cluster_cells) < MIN_CELLS) {

      cat(
        "Skipping cluster:",
        cl,
        "- fewer than",
        MIN_CELLS,
        "cells\n"
      )

      next
    }

    ########################################################
    # BUILD PSEUDOBULKS
    #
    # One column = one biological sample
    # within this cluster
    ########################################################

    cluster_meta <- meta[
      cluster_cells,
      ,
      drop = FALSE
    ]

    sample_ids <- unique(
      cluster_meta[[SAMPLE_COL]]
    )

    pb_list <- list()
    pb_meta_list <- list()

 for (sid in sample_ids) {
      these_cells <- rownames(
        cluster_meta[
          cluster_meta[[SAMPLE_COL]] == sid,
          ,
          drop = FALSE
        ]
      )

      n_cells_available <- length(
        these_cells
      )

      ######################################################
      # SUBSAMPLING
      #
      # FSC / L / S are randomly downsampled to exactly
      # the number of Cluster_G cells in the same sample.
      ######################################################

      if (
        cl %in% SUBSAMPLE_CLUSTERS
      ) {

        if (
          !sid %in% names(reference_cells)
        ) {

          stop(
            "Sample ",
            sid,
            " has no Cluster_G cells. ",
            "Cannot perform matched subsampling for ",
            cl,
            "."
          )
        }

        target_n <- reference_cells[
          sid
        ]

        if (
          n_cells_available < target_n
        ) {

          stop(
            "Cannot subsample ",
            cl,
            " in sample ",
            sid,
            ": ",
            n_cells_available,
            " cells available but Cluster_G has ",
            target_n,
            " cells."
          )
        }

        these_cells <- sample(
          these_cells,
          size = target_n,
          replace = FALSE
        )

        n_cells <- length(
          these_cells
        )

        cat(
          "SUBSAMPLED",
          sid,
          "| cluster:",
          cl,
          "| available:",
          n_cells_available,
          "| target = Cluster_G:",
          target_n,
          "| selected:",
          n_cells,
          "\n"
        )

      } else {

        n_cells <- n_cells_available

        cat(
          "Including",
          sid,
          "| cluster:",
          cl,
          "| cells:",
          n_cells,
          "\n"
        )
      }
      if (n_cells < MIN_CELLS) {

        cat(
          "Excluding",
          sid,
          "from",
          cl,
          ":",
          n_cells,
          "cells <",
          MIN_CELLS,
          "\n"
        )

        next
      }

      stage <- unique(
        cluster_meta[
          these_cells,
          STAGE_COL
        ]
      )


      if (length(stage) != 1) {

        stop(
          "Sample ",
          sid,
          "has multiple stages."
        )
      }

stage <- ifelse(
  stage %in% c("D2", "PM"),
  "D2_PM",
  stage
)

      cat(
        "Including",
        sid,
        "| stage:",
        stage,
        "| cells:",
        n_cells,
        "\n"
      )

      summed <- Matrix::rowSums(
        counts[
          ,
          these_cells,
          drop = FALSE
        ]
      )

      pb_list[[sid]] <- summed

      pb_meta_list[[sid]] <- data.frame(
        sample = sid,
        orig.ident = stage,
        n_cells = n_cells,
        stringsAsFactors = FALSE
      )
    }

    ########################################################
    # CHECK PSEUDOBULK
    ########################################################

    if (length(pb_list) == 0) {

      cat(
        "No valid pseudobulk samples for",
        cl,
        "\n"
      )

      next
    }

    pb_counts <- do.call(
      cbind,
      pb_list
    )

    pb_meta <- do.call(
      rbind,
      pb_meta_list
    )

    rownames(pb_meta) <- pb_meta$sample

    ########################################################
    # REMOVE UNUSED STAGES
    ########################################################

    pb_meta$orig.ident <- factor(
      pb_meta$orig.ident,
      levels = requested.stages
    )

    ########################################################
    # COMPARISONS
    ########################################################

    for (comparison in COMPARISONS) {

      condition_A <- comparison[1]
      condition_B <- comparison[2]

      comparison_name <- paste0(
        condition_A,
        "_vs_",
        condition_B
      )

      cat("\n")
      cat(
        "--------------------------------------------------------\n"
      )

      cat(
        cl,
        ":",
        comparison_name,
        "\n"
      )

      ######################################################
      # SUBSET SAMPLES
      ######################################################

      keep_samples <- rownames(
        pb_meta
      )[
        pb_meta$orig.ident %in%
          c(
            condition_A,
            condition_B
          )
      ]

      sub_meta <- pb_meta[
        keep_samples,
        ,
        drop = FALSE
      ]

      sub_counts <- pb_counts[
        ,
        keep_samples,
        drop = FALSE
      ]

      sub_meta$orig.ident <- droplevels(
        factor(
          sub_meta$orig.ident,
          levels = c(
            condition_B,
            condition_A
          )
        )
      )

      ######################################################
      # REPLICATE CHECK
      ######################################################

      rep_table <- table(
        sub_meta$orig.ident
      )

            ######################################################
      # CHECK CONDITION PRESENCE
      ######################################################

      n_A <- sum(
        sub_meta$orig.ident == condition_A
      )

      n_B <- sum(
        sub_meta$orig.ident == condition_B
      )

      if (
        n_A == 0 ||
        n_B == 0
      ) {

        cat(
          "SKIPPED:",
          comparison_name,
          "- cluster absent from one or both conditions.\n"
        )

        comparison_summary[[
          paste(
            cl,
            comparison_name,
            sep = "__"
          )
        ]] <- data.frame(
          assay = assay_name,
          cluster = cl,
          comparison = comparison_name,
          n_A = n_A,
          n_B = n_B,
          tested_features = 0,
          significant_features = 0,
          stringsAsFactors = FALSE
        )

        next
      }


      cat(
        "Replicates:",
        condition_A,
        "=",
        ifelse(
          condition_A %in% names(rep_table),
          rep_table[condition_A],
          0
        ),
        "|",
        condition_B,
        "=",
        ifelse(
          condition_B %in% names(rep_table),
          rep_table[condition_B],
          0
        ),
        "\n"
      )

      if (
        any(
          rep_table <
            MIN_REPLICATES
        )
      ) {

        cat(
          "SKIPPED:",
          comparison_name,
          "- fewer than",
          MIN_REPLICATES,
          "biological replicates in one condition.\n"
        )

        comparison_summary[[
            paste(
              cl,
              comparison_name,
              sep = "__"
            )
          ]
        ] <- data.frame(
          assay = assay_name,
          cluster = cl,
          comparison = comparison_name,
          n_A = ifelse(
            condition_A %in% names(rep_table),
            rep_table[condition_A],
            0
          ),
          n_B = ifelse(
            condition_B %in% names(rep_table),
            rep_table[condition_B],
            0
          ),
          tested_features = 0,
          significant_features = 0,
          stringsAsFactors = FALSE
        )

        next
      }

      ######################################################
      # LOW COUNT FILTER
      #
      # Keep features with at least 10 counts
      # in at least MIN_REPLICATES samples.
      ######################################################

      keep_features <-
        rowSums(
          sub_counts >= 10
        ) >= MIN_REPLICATES

      sub_counts <- sub_counts[
        keep_features,
        ,
        drop = FALSE
      ]

      cat(
        "Features retained:",
        nrow(sub_counts),
        "\n"
      )

      if (nrow(sub_counts) < 10) {

  cat(
    "SKIPPED:",
    comparison_name,
    "- only",
    nrow(sub_counts),
    "features remain after filtering.\n"
  )

  comparison_summary[[
    paste(
      cl,
      comparison_name,
      sep = "__"
    )
  ]] <- data.frame(
    assay = assay_name,
    cluster = cl,
    comparison = comparison_name,
    n_A = rep_table[condition_A],
    n_B = rep_table[condition_B],
    tested_features = nrow(sub_counts),
    significant_features = 0,
    stringsAsFactors = FALSE
  )

  next
}

      ######################################################
      # DESEQ2
      ######################################################

      dds <- DESeqDataSetFromMatrix(
        countData = round(
        sub_counts
       ),
        colData = sub_meta,
        design = ~ orig.ident
         )


dds <- estimateSizeFactors(
  dds,
  type = "poscounts"
)
      ######################################################
      # DESEQ2
      ######################################################

      dds <- DESeq(
        dds,
        fitType = "local",
        quiet = TRUE
      )

      ######################################################
      # RESULTS
      ######################################################

      res <- results(
        dds,
        contrast = c(
          "orig.ident",
          condition_A,
          condition_B
        ),
        alpha = FDR_THRESHOLD
      )

      res <- as.data.frame(
        res
      )

      res$feature <- rownames(
        res
      )

      ######################################################
      # SIGNIFICANCE
      ######################################################

      res$significant <- (
        !is.na(res$padj) &
        res$padj < FDR_THRESHOLD &
        !is.na(res$log2FoldChange) &
        abs(
          res$log2FoldChange
        ) >= LFC_THRESHOLD
      )

      sig <- res[
        res$significant,
        ,
        drop = FALSE
      ]

      ######################################################
      # SORT
      ######################################################

      res <- res[
        order(
          res$padj,
          na.last = TRUE
        ),
        ,
        drop = FALSE
      ]

      ######################################################
      # SAVE FULL RESULT
      ######################################################

      outfile <- file.path(
        assay_outdir,
        paste0(
          "DESeq2_",
          assay_name,
          "_",
          cl,
          "_",
          comparison_name,
          ".csv"
        )
      )

      write.csv(
        res,
        outfile,
        row.names = FALSE
      )

      ######################################################
      # SUMMARY
      ######################################################

      n_sig <- nrow(
        sig
      )

      cat(
        "SIGNIFICANT:",
        n_sig,
        "\n"
      )

      comparison_summary[[
          paste(
            cl,
            comparison_name,
            sep = "__"
          )
        ]
      ] <- data.frame(
        assay = assay_name,
        cluster = cl,
        comparison = comparison_name,
        n_A = rep_table[condition_A],
        n_B = rep_table[condition_B],
        tested_features = nrow(res),
        significant_features = n_sig,
        stringsAsFactors = FALSE
      )

      ######################################################
      # STORE UNIQUE FEATURES
      ######################################################

      if (n_sig > 0) {

        if (
          is.null(
            unique_features_by_cluster[[cl]]
          )
        ) {

          unique_features_by_cluster[[cl]] <-
            character(0)
        }

        unique_features_by_cluster[[cl]] <-
          unique(
            c(
              unique_features_by_cluster[[cl]],
              sig$feature
            )
          )
      }
    }
  }

  ##########################################################
  # COMPARISON SUMMARY TABLE
  ##########################################################

  summary_df <- if (
    length(comparison_summary) > 0
  ) {

    do.call(
      rbind,
      comparison_summary
    )

  } else {

    data.frame()
  }

  rownames(
    summary_df
  ) <- NULL

  write.csv(
    summary_df,
    file.path(
      assay_outdir,
      paste0(
        "03.5_",
        assay_name,
        "_stage_comparison_summary.csv"
      )
    ),
    row.names = FALSE
  )

  ##########################################################
  # UNIQUE FEATURES ACROSS ALL FOUR COMPARISONS
  ##########################################################

  unique_summary <- data.frame(
    assay = character(),
    cluster = character(),
    total_unique_features = integer(),
    stringsAsFactors = FALSE
  )

  if (
    length(
      unique_features_by_cluster
    ) > 0
  ) {

    for (cl in names(
      unique_features_by_cluster
    )) {

      features <- unique_features_by_cluster[[cl]]

      unique_summary <- rbind(
        unique_summary,
        data.frame(
          assay = assay_name,
          cluster = cl,
          total_unique_features =
            length(
              features
            ),
          stringsAsFactors = FALSE
        )
      )
    }
  }

  write.csv(
    unique_summary,
    file.path(
      assay_outdir,
      paste0(
        "03.5_",
        assay_name,
        "_UNIQUE_total.csv"
      )
    ),
    row.names = FALSE
  )

  ##########################################################
  # CONSOLE TABLE
  ##########################################################

  cat("\n")
  cat("========================================================\n")
  cat(
    "SUMMARY:",
    assay_name,
    "\n"
  )
  cat("========================================================\n")

  if (
    nrow(summary_df) > 0
  ) {

    print(
      summary_df[
        ,
        c(
          "cluster",
          "comparison",
          "n_A",
          "n_B",
          "tested_features",
          "significant_features"
        ),
        drop = FALSE
      ],
      row.names = FALSE
    )
  }

  cat("\n")
  cat("========================================================\n")
  cat(
    "UNIQUE FEATURES ACROSS ALL FOUR COMPARISONS:",
    assay_name,
    "\n"
  )
  cat("========================================================\n")

  if (
    nrow(unique_summary) > 0
  ) {

    print(
      unique_summary[
        ,
        c(
          "cluster",
          "total_unique_features"
        ),
        drop = FALSE
      ],
      row.names = FALSE
    )
  }

  return(
    list(
      summary = summary_df,
      unique = unique_summary
    )
  )
}

############################################################
# RUN RNA
############################################################

cat("\n")
cat("========================================================\n")
cat("STARTING RNA\n")
cat("========================================================\n")

rna_results <- run_stage_analysis(
  obj,
  "SoupXRNA"
)


############################################################
# FINAL COMBINED TABLE
############################################################

cat("\n")
cat("========================================================\n")
cat("FINAL COMBINED SUMMARY\n")
cat("========================================================\n")

final_summary <- rbind(
  rna_results$unique,
  atac_results$unique
)

write.csv(
  final_summary,
  file.path(
    OUTDIR,
    "03.5_FINAL_unique_DEG_by_cluster.csv"
  ),
  row.names = FALSE
)

print(
  final_summary,
  row.names = FALSE
)

cat("\n")
cat("========================================================\n")
cat("03.5 FINISHED\n")
cat("========================================================\n")