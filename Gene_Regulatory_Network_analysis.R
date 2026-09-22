
suppressPackageStartupMessages({
  library(scMEGA)
  library(Seurat)
  library(Signac)
  library(ArchR)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(IRanges)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(purrr)
  library(ComplexHeatmap)
  library(circlize)
  library(Matrix)
  library(igraph)
  library(tidygraph)
  library(ggraph)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(textshape)
  library(progress)
  library(TxDb.Rnorvegicus.UCSC.rn7.refGene)
  library(org.Rn.eg.db)
  library(GenomeInfoDb)
})

set.seed(42)


input_file <- "/shared/projects/femcycle/PG_multiome_cluster/results/objects/PGintegrated73.reclustered.linkpeaks.trajectory.chromvar.rds"

output_dir <- "/shared/projects/femcycle/PG_multiome_cluster/results/GRN"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tf_assay   <- "chromvarPeaks"
rna_assay  <- "SoupXRNA"
atac_assay <- "peaks"

trajectory_name <- "Trajectory"

cells_per_bin <- 10

time_min <- 1
time_max <- 100

tf_p_cutoff   <- 0.05
tf_cor_cutoff <- 0.1

p2g_cor_cutoff <- 0
p2g_fdr_cutoff <- 1e-4

tf_gene_fdr_cutoff <- 1


p2g_distance_cutoff <- 0

genome <- "rn7"

tad_file <- "/shared/projects/femcycle/PG_multiome_cluster/data/TADrn7.bed"

seqlevelsStyle(
    TxDb.Rnorvegicus.UCSC.rn7.refGene
) <- "UCSC"

geneAnnotationRN7 <- createGeneAnnotation(
    TxDb = TxDb.Rnorvegicus.UCSC.rn7.refGene,
    OrgDb = org.Rn.eg.db
)

integrated73 <- readRDS(input_file)

cat("Objet chargé.\n")
cat("Cellules totales :", ncol(integrated73), "\n")
cat(
  "Assays :",
  paste(Seurat::Assays(integrated73), collapse = ", "),
  "\n"
)

if (!"traj_cluster" %in% colnames(integrated73@meta.data)) {
  stop("La colonne 'traj_cluster' est absente de meta.data.")
}

cluster_g_cells <- colnames(integrated73)[
  integrated73$traj_cluster == "Cluster_G"
]

if (length(cluster_g_cells) == 0) {
  stop("Aucune cellule Cluster_G trouvée.")
}

cat(
  "Cellules Cluster_G :",
  length(cluster_g_cells),
  "\n"
)

objG <- subset(
  integrated73,
  cells = cluster_g_cells
)

cat(
  "Objet Cluster_G :",
  ncol(objG),
  "cellules\n"
)

required_assays <- c(
  tf_assay,
  rna_assay,
  atac_assay
)

missing_assays <- setdiff(
  required_assays,
  Seurat::Assays(objG)
)

if (length(missing_assays) > 0) {
  stop(
    "Assays absents de objG : ",
    paste(missing_assays, collapse = ", ")
  )
}

DefaultAssay(objG) <- atac_assay

trajectory_column <- "pt_slingshot2d_Cluster_G"

if (!trajectory_column %in% colnames(objG@meta.data)) {
  stop(
    "La colonne de pseudotemps '",
    trajectory_column,
    "' est absente de objG@meta.data."
  )
}

pseudotime <- objG@meta.data[[trajectory_column]]

names(pseudotime) <- rownames(objG@meta.data)

keep <- is.finite(pseudotime)

if (!any(keep)) {
  stop(
    "Aucune cellule Cluster_G ne possède un pseudotemps valide."
  )
}

objG <- subset(
  objG,
  cells = names(pseudotime)[keep]
)

pseudotime <- pseudotime[colnames(objG)]

cat(
  "Cellules Cluster_G avec pseudotemps valide :",
  length(pseudotime),
  "\n"
)

ordered_cells <- names(
  sort(pseudotime, decreasing = FALSE)
)

objG <- objG[, ordered_cells]

pseudotime <- pseudotime[ordered_cells]

bin_id <- ceiling(
  seq_along(ordered_cells) / cells_per_bin
)

bin_table <- tibble(
  cell = ordered_cells,
  pseudotime = as.numeric(pseudotime),
  bin = bin_id
)

bin_summary <- bin_table %>%
  group_by(bin) %>%
  summarise(
    n_cells = n(),
    pseudotime_min = min(pseudotime),
    pseudotime_max = max(pseudotime),
    pseudotime_mean = mean(pseudotime),
    .groups = "drop"
  )

bin_summary <- bin_summary %>%
  mutate(
    time_point = seq(
      time_min,
      time_max,
      length.out = n()
    )
  )

bin_table <- bin_table %>%
  left_join(
    bin_summary %>%
      dplyr::select(bin, time_point),
    by = "bin"
  )

cell_bin <- setNames(
  bin_table$bin,
  bin_table$cell
)

cell_time_point <- setNames(
  bin_table$time_point,
  bin_table$cell
)

motif_object <- objG[["peaks"]]@motifs

motif_ids <- rownames(
  objG[["chromvarPeaks"]]
)

motif_names <- motif_object@motif.names

if (is.null(names(motif_names))) {
  stop(
    "Les noms des motifs ne permettent pas de construire ",
    "la correspondance JASPAR -> TF."
  )
}

motif_map <- tibble(
  motif_id = names(motif_names),
  tf_label = as.character(motif_names)
) %>%
  filter(
    motif_id %in% motif_ids
  )

if (nrow(motif_map) == 0) {
  stop(
    "Aucune correspondance entre les motifs de peaks ",
    "et chromvarPeaks."
  )
}


motif_map <- motif_map %>%
  mutate(
    tf_label = stringr::str_to_title(tf_label),
    tf_label = stringr::str_replace(
      tf_label,
      "\\.[0-9]+$",
      ""
    )
  )


motif_components <- motif_map %>%
  tidyr::separate_rows(
    tf_label,
    sep = "::"
  ) %>%
  rename(
    tf_component = tf_label
  )

motif_components <- motif_components %>%
  mutate(
    tf_component = stringr::str_trim(tf_component)
  )


cat(
  "Motifs ChromVAR :",
  length(motif_ids),
  "\n"
)

cat(
  "Motifs avec mapping TF :",
  length(unique(motif_components$motif_id)),
  "\n"
)

cat(
  "Motifs sans mapping TF :",
  length(
    setdiff(
      motif_ids,
      unique(motif_components$motif_id)
    )
  ),
  "\n"
)

cat(
  "Nombre de TF distincts associés aux motifs :",
  length(unique(motif_components$tf_component)),
  "\n"
)


DefaultAssay(objG) <- tf_assay

chromvar_matrix <- GetAssayData(
  objG,
  assay = tf_assay,
  layer = "data"
)

chromvar_matrix <- chromvar_matrix[
  intersect(
    motif_ids,
    rownames(chromvar_matrix)
  ),
  ordered_cells,
  drop = FALSE
]

if (nrow(chromvar_matrix) == 0) {
  stop(
    "Aucun motif ChromVAR disponible après filtrage."
  )
}


DefaultAssay(objG) <- rna_assay

rna_matrix <- GetAssayData(
  objG,
  assay = rna_assay,
  layer = "data"
)

rna_matrix <- rna_matrix[
  ,
  ordered_cells,
  drop = FALSE
]


rna_genes <- rownames(rna_matrix)

motif_components <- motif_components %>%
  mutate(
    in_rna = tf_component %in% rna_genes
  )

tf_not_in_rna <- sort(
  unique(
    motif_components$tf_component[
      !motif_components$in_rna
    ]
  )
)

tf_in_rna <- sort(
  unique(
    motif_components$tf_component[
      motif_components$in_rna
    ]
  )
)

cat(
  "TF distincts associés aux motifs :",
  length(
    unique(motif_components$tf_component)
  ),
  "\n"
)

cat(
  "TF présents dans RNA :",
  length(tf_in_rna),
  "\n"
)

cat(
  "TF absents de RNA :",
  length(tf_not_in_rna),
  "\n"
)


aggregate_by_bin <- function(
    matrix,
    bin_vector,
    FUN = mean
) {

  stopifnot(
    length(bin_vector) == ncol(matrix)
  )

  bin_levels <- sort(
    unique(bin_vector)
  )

  out <- vapply(
    bin_levels,
    function(b) {

      cells <- which(
        bin_vector == b
      )

      apply(
        matrix[, cells, drop = FALSE],
        1,
        FUN,
        na.rm = TRUE
      )
    },
    numeric(nrow(matrix))
  )

  out <- as.matrix(out)

  rownames(out) <- rownames(matrix)

  colnames(out) <- paste0(
    "bin_",
    bin_levels
  )

  out
}


bin_vector <- cell_bin[
  colnames(chromvar_matrix)
]


tf_activity_bins <- aggregate_by_bin(
  chromvar_matrix,
  bin_vector
)


rna_bin_vector <- cell_bin[
  colnames(rna_matrix)
]

rna_expression_bins <- aggregate_by_bin(
  rna_matrix,
  rna_bin_vector
)

cat(
  "Dimensions TF activity binée :",
  paste(
    dim(tf_activity_bins),
    collapse = " x "
  ),
  "\n"
)

cat(
  "Dimensions RNA expression binée :",
  paste(
    dim(rna_expression_bins),
    collapse = " x "
  ),
  "\n"
)


tf_monomer_map <- motif_map %>%
  filter(
    !stringr::str_detect(tf_label, "::")
  ) %>%
  filter(
    tf_label %in% rownames(rna_expression_bins)
  )

if (nrow(tf_monomer_map) == 0) {
  stop(
    "Aucun motif monomère ne correspond à un TF présent dans RNA."
  )
}

tf_candidates <- sort(
  unique(tf_monomer_map$tf_label)
)

cat(
  "TF candidats :",
  length(tf_candidates),
  "\n"
)


tf_monomer_map <- tf_monomer_map %>%
  filter(
    motif_id %in% rownames(tf_activity_bins)
  )

tf_activity_candidates <- tf_activity_bins[
  tf_monomer_map$motif_id,
  ,
  drop = FALSE
]

rownames(tf_activity_candidates) <- tf_monomer_map$tf_label


tf_expression_candidates <- rna_expression_bins[
  tf_candidates,
  ,
  drop = FALSE
]


n_bins <- ncol(
  tf_activity_bins
)

cat(
  "Nombre de bins trajectoire :",
  n_bins,
  "\n"
)



correlation_motif_map <- motif_map %>%
  mutate(
    tf_component = strsplit(
      tf_label,
      "::",
      fixed = TRUE
    )
  ) %>%
  tidyr::unnest(
    tf_component
  ) %>%
  mutate(
    tf_component = stringr::str_trim(
      tf_component
    )
  ) %>%
  filter(
    tf_component %in%
      rownames(
        rna_expression_bins
      )
  ) %>%
  dplyr::select(
    motif_id,
    tf_label,
    tf = tf_component
  )


cat(
  "Relations motif x TF soumises à GetCorrelation() :",
  nrow(correlation_motif_map),
  "\n"
)

cat(
  "Motifs représentés :",
  n_distinct(
    correlation_motif_map$motif_id
  ),
  "\n"
)

cat(
  "TF composants représentés :",
  n_distinct(
    correlation_motif_map$tf
  ),
  "\n"
)


correlation_activity <- tf_activity_bins[
  correlation_motif_map$motif_id,
  ,
  drop = FALSE
]



correlation_expression <- rna_expression_bins[
  correlation_motif_map$tf,
  ,
  drop = FALSE
]


rownames(
  correlation_expression
) <- correlation_motif_map$motif_id


if (
  !identical(
    rownames(
      correlation_activity
    ),
    rownames(
      correlation_expression
    )
  )
) {

  stop(
    "Les matrices ChromVAR et RNA ne correspondent pas ",
    "motif par motif."
  )

}


if (
  !identical(
    colnames(
      correlation_activity
    ),
    colnames(
      correlation_expression
    )
  )
) {

  stop(
    "Les bins ChromVAR et RNA ne correspondent pas."
  )

}



trajMM_correlation <-
  SummarizedExperiment::SummarizedExperiment(
    assays = list(
      mat = as.matrix(
        correlation_activity
      )
    )
  )


trajRNA_correlation <-
  SummarizedExperiment::SummarizedExperiment(
    assays = list(
      mat = as.matrix(
        correlation_expression
      )
    )
  )


correlation_raw <- scMEGA::GetCorrelation(
  trajectory1 = trajMM_correlation,
  trajectory2 = trajRNA_correlation
)



tf_correlation_results <- correlation_raw %>%
  tibble::rownames_to_column(
    var = "motif_id"
  ) %>%
  dplyr::select(
    motif_id,
    correlation,
    p_value,
    adj_p
  ) %>%
  rename(
    fdr = adj_p
  ) %>%
  left_join(
    correlation_motif_map %>%
      dplyr::select(
        motif_id,
        tf_label,
        tf
      ),
    by = "motif_id"
  ) %>%
  dplyr::select(
    motif_id,
    tf_label,
    tf,
    correlation,
    p_value,
    fdr
  )


cat(
  "Relations motif x TF avec corrélation valide :",
  sum(
    is.finite(
      tf_correlation_results$correlation
    )
  ),
  "\n"
)

cat(
  "TF distincts avec corrélation valide :",
  tf_correlation_results %>%
    filter(
      is.finite(
        correlation
      )
    ) %>%
    pull(tf) %>%
    unique() %>%
    length(),
  "\n"
)


dimer_correlation_check <-
  tf_correlation_results %>%
  filter(
    str_detect(
      tf_label,
      fixed("::")
    )
  ) %>%
  group_by(
    motif_id,
    tf_label
  ) %>%
  summarise(
    n_components = n(),
    n_significant = sum(
      is.finite(fdr) &
      fdr < 0.05
    ),
    .groups = "drop"
  )


if (
  any(
    dimer_correlation_check$n_components != 2
  )
) {

  stop(
    "Un motif dimérique ne possède pas exactement ",
    "deux composants TF."
  )

}


cat(
  "Motifs dimériques testés :",
  nrow(
    dimer_correlation_check
  ),
  "\n"
)


motif_selection_status <- tf_correlation_results %>%
  group_by(
    motif_id,
    tf_label
  ) %>%
  summarise(
    all_components_pass =
      all(
        is.finite(correlation) &
        is.finite(fdr) &
        correlation > 0.1 &
        fdr < 0.05
      ),
    .groups = "drop"
  )


selected_motif_ids <- motif_selection_status %>%
  filter(
    all_components_pass
  ) %>%
  pull(
    motif_id
  )


selected_tf_motifs_all <- tf_correlation_results %>%
  filter(
    motif_id %in% selected_motif_ids,
    is.finite(correlation),
    is.finite(fdr),
    correlation > 0.1,
    fdr < 0.05
  )


cat(
  "Motifs satisfaisant correlation > 0.1 et FDR < 0.05 :",
  n_distinct(
    selected_tf_motifs_all$motif_id
  ),
  "\n"
)

cat(
  "Relations motif-TF satisfaisant les critères :",
  nrow(
    selected_tf_motifs_all
  ),
  "\n"
)

selected_tf_motifs <- selected_tf_motifs_all %>%
  group_by(
    tf
  ) %>%
  slice_max(
    order_by = correlation,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(
    desc(correlation)
  )


selected_tfs <- sort(
  unique(
    selected_tf_motifs$tf
  )
)


cat(
  "Motifs TF sélectionnés pour les heatmaps :",
  nrow(
    selected_tf_motifs
  ),
  "\n"
)

cat(
  "TF sélectionnés pour les heatmaps :",
  length(
    selected_tfs
  ),
  "\n"
)

heatmap_motif_map <- selected_tf_motifs %>%
  dplyr::select(
    motif_id,
    tf_label,
    tf,
    correlation,
    fdr
  )

cat(
  "Motifs retenus pour les heatmaps (FDR < 0.05) :",
  nrow(heatmap_motif_map),
  "\n"
)

cat(
  "TF distincts représentés :",
  n_distinct(heatmap_motif_map$tf),
  "\n"
)



tf_expression_heatmap_raw <- rna_expression_bins[
  heatmap_motif_map$tf,
  ,
  drop = FALSE
]

rownames(tf_expression_heatmap_raw) <-
  heatmap_motif_map$motif_id


tf_activity_heatmap_raw <- tf_activity_bins[
  heatmap_motif_map$motif_id,
  ,
  drop = FALSE
]


if(
  !identical(
    rownames(tf_activity_heatmap_raw),
    rownames(tf_expression_heatmap_raw)
  )
){

  stop(
    "Les matrices TF activity et TF expression ",
    "ne correspondent pas motif par motif."
  )

}

cat(
  "Correspondance motif activity / expression : OK\n"
)


smooth_trajectory <- function(
    x,
    window = 7
){

  stats::filter(
    x,
    rep(
      1 / window,
      window
    ),
    sides = 2
  ) %>%
    as.numeric()

}



get_expression_peak_time <- function(
    x,
    window = 7,
    top_fraction = 0.95
){

  x <- as.numeric(
    x
  )

  bins <- seq_along(
    x
  )

  smooth_x <- smooth_trajectory(
    x,
    window = window
  )

  valid <- is.finite(
    smooth_x
  )

  if(
    sum(valid) < 3
  ){

    return(
      NA_real_
    )

  }

  x_valid <- smooth_x[
    valid
  ]

  bins_valid <- bins[
    valid
  ]

  xmin <- min(
    x_valid
  )

  xmax <- max(
    x_valid
  )

  amplitude <- xmax - xmin

  if(
    !is.finite(amplitude) ||
    amplitude <= 0
  ){

    return(
      NA_real_
    )

  }

  threshold <-
    xmin +
    top_fraction * amplitude

  top <- smooth_x >= threshold

  top[
    !is.finite(top)
  ] <- FALSE


  r <- rle(
    top
  )

  ends <- cumsum(
    r$lengths
  )

  starts <- c(
    1,
    head(
      ends,
      -1
    ) + 1
  )

  plateau_table <- data.frame(
    start = starts[r$values],
    end = ends[r$values]
  )

  if(
    nrow(plateau_table) == 0
  ){

    return(
      bins[
        which.max(
          smooth_x
        )
      ]
    )

  }



  max_bin <- bins[
    which.max(
      smooth_x
    )
  ]

  plateau <- plateau_table[
    plateau_table$start <= max_bin &
      plateau_table$end >= max_bin,
    ,
    drop = FALSE
  ]

  if(
    nrow(plateau) == 0
  ){

    return(
      max_bin
    )

  }

  start <- plateau$start[1]
  end <- plateau$end[1]

  idx <- start:end

  weights <- smooth_x[
    idx
  ]

  weights <- weights -
    min(
      weights,
      na.rm = TRUE
    )

  if(
    sum(weights, na.rm = TRUE) <= 0
  ){

    return(
      mean(idx)
    )

  }

  weighted.mean(
    idx,
    weights,
    na.rm = TRUE
  )

}


smooth_window_bins <- max(
  3,
  round(
    0.07 * ncol(
      tf_expression_heatmap_raw
    )
  )
)

cat(
  "Fenêtre de lissage :",
  smooth_window_bins,
  "bins (7% de la trajectoire)\n"
)


tf_expression_peak <- apply(
  tf_expression_heatmap_raw,
  1,
  get_expression_peak_time,
  window = smooth_window_bins,
  top_fraction = 0.95
)


if(
  any(
    !is.finite(
      tf_expression_peak
    )
  )
){

  stop(
    "Certains motifs n'ont pas de position temporelle ",
    "d'expression définissable."
  )

}


motif_order <- order(
  tf_expression_peak,
  heatmap_motif_map$tf,
  heatmap_motif_map$motif_id
)

heatmap_motif_map <- heatmap_motif_map[
  motif_order,
  ,
  drop = FALSE
]

tf_expression_peak <- tf_expression_peak[
  motif_order
]


tf_expression_heatmap_raw <-
  tf_expression_heatmap_raw[
    heatmap_motif_map$motif_id,
    ,
    drop = FALSE
  ]

tf_activity_heatmap_raw <-
  tf_activity_heatmap_raw[
    heatmap_motif_map$motif_id,
    ,
    drop = FALSE
  ]



make_scaled <- function(
    m
){

  m <- as.matrix(
    m
  )

  if(
    nrow(m) == 0
  ){

    return(
      m
    )

  }

  z <- t(
    scale(
      t(m)
    )
  )

  z

}


tf_expression_heatmap <- make_scaled(
  tf_expression_heatmap_raw
)

tf_activity_heatmap <- make_scaled(
  tf_activity_heatmap_raw
)


tf_expression_heatmap[
  tf_expression_heatmap < -2
] <- -2

tf_expression_heatmap[
  tf_expression_heatmap > 2
] <- 2

tf_activity_heatmap[
  tf_activity_heatmap < -2
] <- -2

tf_activity_heatmap[
  tf_activity_heatmap > 2
] <- 2


if(
  !identical(
    rownames(tf_expression_heatmap),
    rownames(tf_activity_heatmap)
  )
){

  stop(
    "L'ordre des motifs n'est pas identique entre ",
    "expression et activité ChromVAR."
  )

}

cat(
  "Ordre RNA -> ChromVAR : OK\n"
)


tf_time_points <- heatmap_motif_map %>%
  mutate(
    peak_bin = as.numeric(
      tf_expression_peak
    ),
    time_point = (
      peak_bin - 1
    ) /
    (
      ncol(tf_expression_heatmap) - 1
    ) *
    100
  ) %>%
  dplyr::select(
    motif_id,
    tf,
    correlation,
    fdr,
    peak_bin,
    time_point
  )


write.csv(
  tf_activity_heatmap,
  file.path(
    output_dir,
    "TF_activity_heatmap_matrix.csv"
  ),
  quote = FALSE
)

write.csv(
  tf_expression_heatmap,
  file.path(
    output_dir,
    "TF_expression_heatmap_matrix.csv"
  ),
  quote = FALSE
)

write.csv(
  tf_time_points,
  file.path(
    output_dir,
    "TF_heatmap_time_points.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


ExportTFHeatmaps <- function(
    activity_matrix,
    expression_matrix,
    output_dir
){

  if(
    !identical(
      rownames(activity_matrix),
      rownames(expression_matrix)
    )
  ){

    stop(
      "Les deux heatmaps doivent avoir ",
      "exactement le même ordre de motifs."
    )

  }


  col_tf <- circlize::colorRamp2(
    c(
      -2,
      0,
      2
    ),
    c(
      "#2166AC",
      "#F7F7F7",
      "#B2182B"
    )
  )

  col_rna <- circlize::colorRamp2(
    c(
      -2,
      -1,
      0,
      1,
      2
    ),
    c(
      "#0C0786",
      "#5B02A3",
      "#B5367A",
      "#E56B5D",
      "#F0F921"
    )
  )


  ht_tf <- ComplexHeatmap::Heatmap(
    activity_matrix,
    name = "TF activity",
    col = col_tf,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    use_raster = TRUE
  )


  label_width <-
    ComplexHeatmap::max_text_width(
      heatmap_motif_map$tf,
      gp = grid::gpar(
        fontsize = 7
      )
    ) +
    grid::unit(
      3,
      "mm"
    )


  lab_anno <- ComplexHeatmap::rowAnnotation(
    TF = ComplexHeatmap::anno_text(
      heatmap_motif_map$tf,
      which = "row",
      gp = grid::gpar(
        fontsize = 7
      ),
      just = "right",
      location = 0.5,
      rot = 0
    ),
    width = label_width
  )


  ht_rna <- ComplexHeatmap::Heatmap(
    expression_matrix,
    name = "TF expression",
    col = col_rna,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    use_raster = TRUE
  )


  ht <- ht_tf +
    lab_anno +
    ht_rna


  pdf(
    file.path(
      output_dir,
      "TF_activity_expression_heatmaps.pdf"
    ),
    width = 12,
    height = max(
      8,
      0.20 * nrow(activity_matrix)
    )
  )

  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "right"
  )

  dev.off()


  invisible(
    ht
  )

}


tf_heatmap_object <- ExportTFHeatmaps(
  activity_matrix =
    tf_activity_heatmap,
  expression_matrix =
    tf_expression_heatmap,
  output_dir =
    output_dir
)


write.csv(
  heatmap_motif_map,
  file.path(
    output_dir,
    "TF_heatmap_selected_motifs_FDR005.csv"
  ),
  row.names = FALSE,
  quote = FALSE
)


DefaultAssay(objG) <- atac_assay

atac_matrix <- GetAssayData(
    objG,
    assay = atac_assay,
    layer = "data"
)

atac_matrix <- atac_matrix[
    ,
    ordered_cells,
    drop = FALSE
]

if (
    !identical(
        colnames(atac_matrix),
        ordered_cells
    )
) {

    stop(
        "L'ordre des cellules de la matrice ATAC ",
        "ne correspond pas à ordered_cells."
    )

}

atac_bin_vector <- cell_bin[
    colnames(atac_matrix)
]


if (
    !identical(
        names(atac_bin_vector),
        colnames(atac_matrix)
    )
) {

    stop(
        "Le mapping cellule -> bin ATAC est incorrect."
    )

}

atac_expression_bins <- aggregate_by_bin(
    atac_matrix,
    atac_bin_vector
)


cat(
    "Dimensions ATAC binée :",
    paste(
        dim(atac_expression_bins),
        collapse = " x "
    ),
    "\n"
)


if (
    ncol(atac_expression_bins) !=
    ncol(rna_expression_bins)
) {

    stop(
        "Le nombre de bins ATAC (",
        ncol(atac_expression_bins),
        ") ne correspond pas au nombre de bins RNA (",
        ncol(rna_expression_bins),
        ")."
    )

}


if (
    !identical(
        colnames(atac_expression_bins),
        colnames(rna_expression_bins)
    )
) {

    stop(
        "Les noms des bins ATAC et RNA ne correspondent pas."
    )

}


cat(
    "Correspondance des bins ATAC / RNA : OK\n"
)


smooth_peak_gene_window <- 18
smooth_atac_window <- smooth_peak_gene_window

smoothTrajectoryMatrix <- function(
    mat,
    window = 18
) {

    mat <- as.matrix(mat)

    n_bins <- ncol(mat)

    if (window < 1) {
        stop("window doit être >= 1.")
    }

    if (window > n_bins) {
        stop(
            "La fenêtre de lissage (",
            window,
            ") est supérieure au nombre de bins (",
            n_bins,
            ")."
        )
    }

    out <- matrix(
        NA_real_,
        nrow = nrow(mat),
        ncol = ncol(mat),
        dimnames = dimnames(mat)
    )


    half_window <- floor(
        window / 2
    )

    for (i in seq_len(n_bins)) {

        start_bin <- max(
            1,
            i - half_window
        )

        end_bin <- min(
            n_bins,
            i + half_window
        )

        out[, i] <- rowMeans(
            mat[
                ,
                start_bin:end_bin,
                drop = FALSE
            ],
            na.rm = TRUE
        )
    }

    out
}


rna_expression_bins_smoothed <-
    smoothTrajectoryMatrix(
        rna_expression_bins,
        window = smooth_peak_gene_window
    )


cat(
    "Dimensions RNA lissée :",
    paste(
        dim(
            rna_expression_bins_smoothed
        ),
        collapse = " x "
    ),
    "\n"
)

if (
    anyNA(
        rna_expression_bins_smoothed
    )
) {

    stop(
        "Le lissage RNA contient des NA."
    )

}

cat(
    "RNA binée + lissée : OK\n"
)


atac_expression_bins_smoothed <- smoothTrajectoryMatrix(
    atac_expression_bins,
    window = smooth_atac_window
)

cat(
    "Dimensions ATAC lissée :",
    paste(
        dim(
            atac_expression_bins_smoothed
        ),
        collapse = " x "
    ),
    "\n"
)

if (
    anyNA(
        atac_expression_bins_smoothed
    )
) {

    stop(
        "Le lissage ATAC contient des NA."
    )

}

cat(
    "ATAC binée + lissée : OK\n"
)


if (
    !identical(
        colnames(
            rna_expression_bins_smoothed
        ),
        colnames(
            atac_expression_bins_smoothed
        )
    )
) {

    stop(
        "Les bins RNA et ATAC lissés ne correspondent pas."
    )

}

cat(
    "Correspondance RNA / ATAC après lissage : OK\n"
)

cat(
    "Fenêtre de lissage Peak-to-Gene :",
    smooth_peak_gene_window,
    "bins\n"
)

if (
    !identical(
        rownames(atac_expression_bins),
        rownames(atac_expression_bins_smoothed)
    )
) {

    stop(
        "Les peaks ne correspondent plus après lissage."
    )

}


if (
    !identical(
        colnames(atac_expression_bins),
        colnames(atac_expression_bins_smoothed)
    )
) {

    stop(
        "Les bins ATAC ne correspondent plus après lissage."
    )

}


cat(
    "ATAC binée + lissée : OK\n"
)

cat(
    "Fenêtre de lissage ATAC :",
    smooth_atac_window,
    "bins\n"
)


p2g_max_dist <- 50000000

p2g_cor_cutoff <- 0.7

p2g_fdr_cutoff <- 1e-4

tad_file <- "/shared/projects/femcycle/PG_multiome_cluster/data/TADrn7.bed"


if (!file.exists(tad_file)) {

    stop(
        "Fichier TAD introuvable : ",
        tad_file,
        "\n",
        "Transférer TADrn7.bed dans ",
        file.path(project_dir, "data")
    )

}

cat(
    "TAD utilisé : ",
    tad_file,
    "\n"
)

rna_peak_gene_bins <-  rna_expression_bins_smoothed

atac_peak_gene_bins <- atac_expression_bins_smoothed

cat(
    "RNA utilisé pour Peak-to-Gene : biné + lissé\n"
)

cat(
    "ATAC utilisé pour Peak-to-Gene : biné + lissé\n"
)


cat(
    "Dimensions ATAC binée + lissée : ",
    paste(
        dim(atac_peak_gene_bins),
        collapse = " x "
    ),
    "\n"
)

cat(
    "Dimensions RNA binée : ",
    paste(
        dim(rna_peak_gene_bins),
        collapse = " x "
    ),
    "\n"
)

if (
    !identical(
        colnames(atac_peak_gene_bins),
        colnames(rna_peak_gene_bins)
    )
) {

    stop(
        "Les bins ATAC et RNA ne correspondent pas."
    )

}

cat(
    "Correspondance ATAC / RNA : OK\n"
)



PeakToGene <- function(
    peak.mat,
    gene.mat,
    genome = "rn7",
    max.dist = 50000000,
    method = "correlation",
    tad.file = NULL,
    workers = 1
) {

    if (method != "correlation") {
        stop(
            "Seule la méthode 'correlation' est supportée."
        )
    }

    if (genome != "rn7") {
        stop(
            "Cette version de PeakToGene est prévue pour rn7."
        )
    }

    if (!exists("geneAnnotationRN7")) {
        stop(
            "geneAnnotationRN7 n'existe pas dans la session."
        )
    }



    genes <- geneAnnotationRN7$genes

    gene_symbols <- as.character(
        genes$symbol
    )

    keep_gene <- (
        gene_symbols %in%
        rownames(gene.mat)
    )

    genes <- genes[keep_gene]

    gene_symbols <- gene_symbols[
        keep_gene
    ]

    gene_mat <- gene.mat[
        gene_symbols,
        ,
        drop = FALSE
    ]

    # Une seule occurrence par symbole
    keep_unique <- !duplicated(
        gene_symbols
    )

    genes <- genes[
        keep_unique
    ]

    gene_symbols <- gene_symbols[
        keep_unique
    ]

    gene_mat <- gene_mat[
        keep_unique,
        ,
        drop = FALSE
    ]



    gene_start <- ifelse(
        as.character(
            strand(genes)
        ) == "+",
        start(genes),
        end(genes)
    )

    gene_gr <- GRanges(
        seqnames = seqnames(genes),
        ranges = IRanges(
            start = gene_start,
            width = 1
        ),
        strand = strand(genes),
        gene = gene_symbols
    )


    peak_names <- rownames(
        peak.mat
    )

    peak_split <- stringr::str_split_fixed(
        peak_names,
        "-",
        3
    )

    peak_gr <- GRanges(
        seqnames = peak_split[, 1],
        ranges = IRanges(
            start = as.numeric(
                peak_split[, 2]
            ),
            end = as.numeric(
                peak_split[, 3]
            )
        )
    )

    if (!is.null(tad.file)) {

        if (!file.exists(tad.file)) {
            stop(
                "Fichier TAD introuvable : ",
                tad.file
            )
        }

        tads <- read.table(
            tad.file,
            header = FALSE,
            sep = "\t",
            stringsAsFactors = FALSE
        )

        tads <- tads[, 1:3]

        colnames(tads) <- c(
            "chr",
            "start",
            "end"
        )

        tads_gr <- makeGRangesFromDataFrame(
            tads,
            seqnames.field = "chr",
            start.field = "start",
            end.field = "end"
        )

    } else {

        tads_gr <- NULL

    }

    gene_tad <- rep(
        NA_integer_,
        length(gene_gr)
    )

    peak_tad <- rep(
        NA_integer_,
        length(peak_gr)
    )

    if (!is.null(tads_gr)) {

        gene_hits <- findOverlaps(
            gene_gr,
            tads_gr,
            ignore.strand = TRUE
        )

        peak_hits <- findOverlaps(
            peak_gr,
            tads_gr,
            ignore.strand = TRUE
        )

        gene_tad[
            queryHits(gene_hits)
        ] <- subjectHits(
            gene_hits
        )

        peak_tad[
            queryHits(peak_hits)
        ] <- subjectHits(
            peak_hits
        )
    }

    cat(
        "Génération des couples peak-gene candidats...\n"
    )

    gene_windows <- resize(
        gene_gr,
        width = 2 * max.dist + 1,
        fix = "center"
    )

    overlaps <- findOverlaps(
        gene_windows,
        peak_gr,
        ignore.strand = TRUE
    )

    candidate <- data.frame(
        gene_idx = queryHits(
            overlaps
        ),
        peak_idx = subjectHits(
            overlaps
        )
    )

    if (nrow(candidate) == 0) {
        stop(
            "Aucun couple peak-gene trouvé."
        )
    }
    candidate$distance <- distance(
        gene_gr[
            candidate$gene_idx
        ],
        peak_gr[
            candidate$peak_idx
        ]
    )

    candidate <- candidate[
        candidate$distance <= max.dist,
        ,
        drop = FALSE
    ]


    candidate$gene_tad <- gene_tad[
        candidate$gene_idx
    ]

    candidate$peak_tad <- peak_tad[
        candidate$peak_idx
    ]

    if (!is.null(tads_gr)) {

        candidate <- candidate[
            !is.na(candidate$gene_tad) &
            !is.na(candidate$peak_tad) &
            candidate$gene_tad ==
            candidate$peak_tad,
            ,
            drop = FALSE
        ]
    }

    cat(
        "Couples après filtrage génomique :",
        format(
            nrow(candidate),
            big.mark = ","
        ),
        "\n"
    )

    if (nrow(candidate) == 0) {
        stop(
            "Aucun couple peak-gene après filtrage TAD/distance."
        )
    }

    atac <- as.matrix(
        peak.mat[
            candidate$peak_idx,
            ,
            drop = FALSE
        ]
    )

    rna <- as.matrix(
        gene_mat[
            candidate$gene_idx,
            ,
            drop = FALSE
        ]
    )


cat(
    "Vérification des trajectoires ATAC/RNA...\n"
)


valid_pair_count <- rowSums(
    is.finite(atac) &
    is.finite(rna)
)

valid <- valid_pair_count >= 3

cat(
    "Couples avec >= 3 bins conjointement finis :",
    sum(valid),
    "\n"
)

cat(
    "Couples insuffisants :",
    sum(!valid),
    "\n"
)

candidate <- candidate[
    valid,
    ,
    drop = FALSE
]

atac <- atac[
    valid,
    ,
    drop = FALSE
]

rna <- rna[
    valid,
    ,
    drop = FALSE
]

valid_pair_count <- valid_pair_count[
    valid
]

if (nrow(candidate) == 0) {

    stop(
        "Aucun couple peak-gene ne possède suffisamment ",
        "de bins ATAC/RNA conjointement finis."
    )

}

    workers <- as.integer(
        workers
    )

    if (
        !is.finite(workers) ||
        workers < 1
    ) {
        workers <- 1
    }

    workers <- min(
        workers,
        parallel::detectCores()
    )

    cat(
        "Calcul des corrélations avec",
        workers,
        "CPU...\n"
    )

correlation_block <- function(
    idx
) {

    x <- atac[
        idx,
        ,
        drop = FALSE
    ]

    y <- rna[
        idx,
        ,
        drop = FALSE
    ]

    cor_values <- numeric(
        nrow(x)
    )

    for (j in seq_len(nrow(x))) {

        ok <- (
            is.finite(x[j, ]) &
            is.finite(y[j, ])
        )

        if (sum(ok) < 3) {

            cor_values[j] <- NA_real_

        } else {

            cor_values[j] <- cor(
                x[j, ok],
                y[j, ok],
                method = "pearson"
            )
        }
    }

    cor_values
}

    # Découpage en blocs
    block_id <- cut(
        seq_len(
            nrow(candidate)
        ),
        breaks = workers,
        labels = FALSE
    )

    index_blocks <- split(
        seq_len(
            nrow(candidate)
        ),
        block_id
    )


    if (workers > 1) {

        cor_blocks <- parallel::mclapply(
            index_blocks,
            correlation_block,
            mc.cores = workers
        )

    } else {

        cor_blocks <- lapply(
            index_blocks,
            correlation_block
        )
    }


    candidate$Correlation <- unlist(
        cor_blocks,
        use.names = FALSE
    )

    n_bins <- ncol(
        atac
    )

    candidate$TStat <- (
        candidate$Correlation /
        sqrt(
            pmax(
                1 -
                candidate$Correlation^2,
                .Machine$double.eps
            ) /
            (
                n_bins - 2
            )
        )
    )

    candidate$Pval <- 2 * pt(
        -abs(
            candidate$TStat
        ),
        df = n_bins - 2
    )

    candidate$FDR <- p.adjust(
        candidate$Pval,
        method = "fdr"
    )

    candidate$gene <- gene_symbols[
        candidate$gene_idx
    ]

    candidate$peak <- peak_names[
        candidate$peak_idx
    ]

    candidate <- candidate[
        is.finite(
            candidate$FDR
        ),
        ,
        drop = FALSE
    ]

    candidate <- candidate[
        order(
            candidate$FDR
        ),
        ,
        drop = FALSE
    ]

    cat(
        "Peak-to-Gene terminé :",
        format(
            nrow(candidate),
            big.mark = ","
        ),
        "relations valides.\n"
    )

    return(
        candidate
    )
}

cat(
    "Calcul des relations Peak-to-Gene...\n"
)

p2g_raw <- PeakToGene(
    peak.mat = atac_peak_gene_bins,
    gene.mat = rna_peak_gene_bins,
    genome = "rn7",
    max.dist = p2g_max_dist,
    method = "correlation",
    tad.file = tad_file,
    workers = as.integer(
        Sys.getenv(
            "SLURM_CPUS_PER_TASK",
            "1"
        )
    )
)

cat(
    "Relations Peak-to-Gene brutes :",
    nrow(p2g_raw),
    "\n"
)

df.p2g <- p2g_raw %>%
    dplyr::filter(
        distance > 0,
        Correlation > p2g_cor_cutoff,
        FDR < p2g_fdr_cutoff
    )


cat(
    "Relations Peak-to-Gene après filtrage :",
    nrow(df.p2g),
    "\n"
)

var_cutoff_gene <- 0.80

n_trajectory_bins <- ncol(
    rna_peak_gene_bins
)

if (
    n_trajectory_bins !=
    ncol(atac_peak_gene_bins)
) {

    stop(
        "RNA et ATAC n'ont pas le même nombre de bins."
    )

}

peak_time_window <- smooth_peak_gene_window


message(
    "Fenêtre utilisée pour les positions temporelles : ",
    peak_time_window,
    " bins"
)

message(
    "Calcul de la position temporelle des gènes RNA ",
    "par sommet / plateau..."
)

genes_needed <- unique(
    df.p2g$gene
)

genes_needed <- intersect(
    genes_needed,
    rownames(
        rna_peak_gene_bins
    )
)

message(
    "Gènes RNA nécessaires pour Peak-to-Gene : ",
    length(genes_needed)
)


if (
    length(genes_needed) == 0
) {

    stop(
        "Aucun gène de df.p2g n'est présent dans ",
        "rna_peak_gene_bins."
    )

}
gene_mat <- rna_peak_gene_bins[
    genes_needed,
    ,
    drop = FALSE
]
gene_has_na <- rowSums(
    is.na(gene_mat)
) > 0


if (
    any(gene_has_na)
) {

    message(
        "Gènes nécessaires contenant des NA supprimés : ",
        sum(gene_has_na)
    )

    gene_mat <- gene_mat[
        !gene_has_na,
        ,
        drop = FALSE
    ]

}

gene_sd <- matrixStats::rowSds(
    gene_mat
)

gene_variable <- is.finite(
    gene_sd
) & gene_sd > 0


message(
    "Gènes nécessaires avec variance > 0 : ",
    sum(gene_variable)
)

message(
    "Gènes nécessaires entièrement nuls supprimés : ",
    sum(!gene_variable)
)


gene_mat <- gene_mat[
    gene_variable,
    ,
    drop = FALSE
]


if (
    nrow(gene_mat) == 0
) {

    stop(
        "Aucun gène nécessaire au Peak-to-Gene ",
        "ne possède une variance non nulle."
    )

}

gene_peak_bin <- apply(
    gene_mat,
    1,
    get_expression_peak_time,
    window = peak_time_window,
    top_fraction = 0.95
)


if (
    any(
        !is.finite(
            gene_peak_bin
        )
    )
) {

    stop(
        "Certains gènes RNA nécessaires au Peak-to-Gene ",
        "n'ont pas de position temporelle valide."
    )

}

gene_time_point <- (
    gene_peak_bin - 1
) /
(
    n_trajectory_bins - 1
) *
100


df_gene_time_point <- data.frame(

    gene =
        names(
            gene_peak_bin
        ),

    peak_bin =
        as.numeric(
            gene_peak_bin
        ),

    gene_time_point =
        as.numeric(
            gene_time_point
        ),

    stringsAsFactors = FALSE
)


message(
    "Temps gènes RNA calculés : ",
    nrow(df_gene_time_point)
)

message(
    "Calcul de la position temporelle des peaks ATAC par sommet / plateau..."
)

peaks_needed <- unique(df.p2g$peak)

peaks_needed <- intersect(
    peaks_needed,
    rownames(atac_peak_gene_bins)
)

message(
    "Peaks ATAC nécessaires pour Peak-to-Gene : ",
    length(peaks_needed)
)

atac_peak_valid <- atac_peak_gene_bins[
    peaks_needed,
    ,
    drop = FALSE
]

peak_sd <- matrixStats::rowSds(
    atac_peak_valid
)

peak_variable <- is.finite(peak_sd) & peak_sd != 0

atac_peak_valid <- atac_peak_valid[
    peak_variable,
    ,
    drop = FALSE
]

message(
    "Peaks ATAC nécessaires avec variance > 0 : ",
    nrow(atac_peak_valid)
)

atac_peak_bin <- apply(
    atac_peak_valid,
    1,
    get_expression_peak_time,
    window = peak_time_window,
    top_fraction = 0.95
)

if (
    any(
        !is.finite(atac_peak_bin)
    )
) {
    stop(
        "Certains peaks ATAC nécessaires au Peak-to-Gene ",
        "n'ont pas de position temporelle valide."
    )
}

atac_time_point <- (
    atac_peak_bin - 1
) /
(
    n_trajectory_bins - 1
) *
100

df_atac_time_point <- data.frame(
    peak = names(atac_peak_bin),
    peak_bin = as.numeric(atac_peak_bin),
    atac_time_point = as.numeric(atac_time_point),
    stringsAsFactors = FALSE
)

if (
    anyDuplicated(
        df_atac_time_point$peak
    ) > 0
) {

    stop(
        "Peaks ATAC dupliqués dans df_atac_time_point."
    )

}


if (
    anyDuplicated(
        df_gene_time_point$gene
    ) > 0
) {

    stop(
        "Gènes dupliqués dans df_gene_time_point."
    )

}


if (
    any(
        !is.finite(
            df_atac_time_point$atac_time_point
        )
    )
) {

    stop(
        "Certains peaks ATAC ont un time point invalide."
    )

}


if (
    any(
        !is.finite(
            df_gene_time_point$gene_time_point
        )
    )
) {

    stop(
        "Certains gènes RNA ont un time point invalide."
    )

}


message(
    "Temps peaks ATAC calculés : ",
    nrow(df_atac_time_point)
)

message(
    "Temps gènes RNA calculés : ",
    nrow(df_gene_time_point)
)

message(
    "Plage temporelle ATAC : ",
    round(
        min(
            df_atac_time_point$atac_time_point
        ),
        2
    ),
    " - ",
    round(
        max(
            df_atac_time_point$atac_time_point
        ),
        2
    )
)

message(
    "Plage temporelle RNA : ",
    round(
        min(
            df_gene_time_point$gene_time_point
        ),
        2
    ),
    " - ",
    round(
        max(
            df_gene_time_point$gene_time_point
        ),
        2
    )
)

message(
    "Exemples des premiers gènes selon leur temps :"
)

print(
    head(
        df_gene_time_point[
            order(
                df_gene_time_point$gene_time_point
            ),
        ],
        10
    )
)

message(
    "Exemples des derniers gènes selon leur temps :"
)

print(
    tail(
        df_gene_time_point[
            order(
                df_gene_time_point$gene_time_point
            ),
        ],
        10
    )
)

df.p2g <- df.p2g %>%

    dplyr::left_join(
        df_atac_time_point,
        by = "peak"
    ) %>%

    dplyr::left_join(
        df_gene_time_point,
        by = "gene"
    ) %>%

    dplyr::mutate(
        time_diff =
            abs(
                atac_time_point -
                gene_time_point
            )
    )


if (
    any(
        !is.finite(
            df.p2g$atac_time_point
        )
    )
) {

    stop(
        "Certaines relations Peak-to-Gene n'ont pas ",
        "de temps ATAC correspondant."
    )

}


if (
    any(
        !is.finite(
            df.p2g$gene_time_point
        )
    )
) {

    stop(
        "Certaines relations Peak-to-Gene n'ont pas ",
        "de temps RNA correspondant."
    )

}


if (
    any(
        !is.finite(
            df.p2g$time_diff
        )
    )
) {

    stop(
        "Certaines relations Peak-to-Gene ont un ",
        "time_diff invalide."
    )

}


message(
    "Correspondance temporelle Peak-to-Gene : OK"
)

message(
    "Relations Peak-to-Gene avec temps : ",
    nrow(df.p2g)
)

message(
    "Distribution de time_diff :"
)

print(
    summary(
        df.p2g$time_diff
    )
)
tf.use <- unique(
    selected_tf_motifs$tf
)

# Vérification dans la matrice RNA réelle
tf.use <- intersect(
    tf.use,
    rownames(rna_expression_bins_smoothed)
)

message(
    "TF utilisés pour TF-Gene : ",
    length(tf.use)
)

if (length(tf.use) == 0) {

    stop(
        "Aucun TF sélectionné n'est présent dans ",
        "rna_expression_bins_smoothed."
    )
}

gene.use <- unique(
    df.p2g$gene
)

gene.use <- intersect(
    gene.use,
    rownames(rna_expression_bins_smoothed)
)

message(
    "Gènes utilisés pour TF-Gene : ",
    length(gene.use)
)


if (length(tf.use) == 0) {

    stop(
        "Aucun TF sélectionné n'est présent dans la matrice TF expression."
    )
}


if (length(gene.use) == 0) {

    stop(
        "Aucun gène Peak-to-Gene n'est présent dans la matrice RNA."
    )
}

tf_motif_map_grn <- selected_tf_motifs_all %>%
    dplyr::select(
        motif_id,
        tf
    ) %>%
    dplyr::distinct()


message(
    "Relations motif-TF utilisées pour TF-Gene : ",
    nrow(tf_motif_map_grn)
)

message(
    "TF distincts utilisés : ",
    dplyr::n_distinct(
        tf_motif_map_grn$tf
    )
)

message(
    "Motifs distincts utilisés : ",
    dplyr::n_distinct(
        tf_motif_map_grn$motif_id
    )
)


tf.use <- intersect(
    unique(
        tf_motif_map_grn$tf
    ),
    rownames(
        rna_expression_bins_smoothed
    )
)


message(
    "TF utilisés pour TF-Gene : ",
    length(tf.use)
)


if (
    length(tf.use) == 0
) {

    stop(
        "Aucun TF sélectionné n'est présent dans ",
        "rna_expression_bins_smoothed."
    )
}

gene.use <- intersect(
    unique(
        df.p2g$gene
    ),
    rownames(
        rna_expression_bins_smoothed
    )
)


message(
    "Gènes utilisés pour TF-Gene : ",
    length(gene.use)
)


if (
    length(gene.use) == 0
) {

    stop(
        "Aucun gène Peak-to-Gene n'est présent dans ",
        "rna_expression_bins_smoothed."
    )
}

gene_expression_for_cor <- rna_expression_bins_smoothed[
    gene.use,
    ,
    drop = FALSE
]


tf_gene_by_motif <- list()


for (
    i in seq_len(
        nrow(tf_motif_map_grn)
    )
) {

    motif_id_i <-
        tf_motif_map_grn$motif_id[i]

    tf_i <-
        tf_motif_map_grn$tf[i]


    if (
        !motif_id_i %in%
        rownames(tf_activity_bins)
    ) {
        next
    }


    motif_activity_i <-
        tf_activity_bins[
            motif_id_i,
            ,
            drop = FALSE
        ]


    cor_i <- cor(
        as.numeric(
            motif_activity_i[1, ]
        ),
        t(
            gene_expression_for_cor
        ),
        use = "pairwise.complete.obs",
        method = "pearson"
    )


    cor_i <- as.numeric(
        cor_i
    )


    names(cor_i) <-
        colnames(
            t(
                gene_expression_for_cor
            )
        )


    tmp <- data.frame(
        tf = tf_i,
        motif_id = motif_id_i,
        gene = gene.use,
        correlation = cor_i,
        stringsAsFactors = FALSE
    )


    tf_gene_by_motif[[i]] <- tmp
}


tf.gene.by.motif <- dplyr::bind_rows(
    tf_gene_by_motif
)

n_cor_bins <- ncol(
    tf_activity_bins
)


tf.gene.by.motif <- tf.gene.by.motif %>%
    dplyr::mutate(
        t_stat =
            correlation /
            sqrt(
                pmax(
                    1 - correlation^2,
                    1e-15
                ) /
                (n_cor_bins - 2)
            ),
        p_value =
            2 *
            pt(
                -abs(t_stat),
                df = n_cor_bins - 2
            )
    ) %>%
    dplyr::group_by(
        tf
    ) %>%
    dplyr::mutate(
        fdr = p.adjust(
            p_value,
            method = "fdr"
        )
    ) %>%
    dplyr::ungroup()

tf.gene.cor <- tf.gene.by.motif %>%
    dplyr::filter(
        is.finite(correlation),
        is.finite(p_value),
        is.finite(fdr)
    ) %>%
    dplyr::group_by(
        tf,
        gene
    ) %>%
    dplyr::slice_max(
        order_by = abs(correlation),
        n = 1,
        with_ties = FALSE
    ) %>%
    dplyr::ungroup()

tf.gene.heatmap <- tf.gene.by.motif %>%
    dplyr::filter(
        is.finite(correlation),
        is.finite(p_value),
        is.finite(fdr)
    ) %>%
    dplyr::group_by(
        tf,
        gene
    ) %>%
    dplyr::filter(
        any(
            fdr < tf_gene_fdr_cutoff
        )
    ) %>%
    dplyr::slice_max(
        order_by = abs(correlation),
        n = 1,
        with_ties = FALSE
    ) %>%
    dplyr::ungroup()


message(
    "Relations TF-Gene retenues pour la heatmap : ",
    nrow(tf.gene.heatmap)
)

message(
    "TF distincts dans la heatmap : ",
    dplyr::n_distinct(
        tf.gene.heatmap$tf
    )
)

message(
    "Gènes distincts dans la heatmap : ",
    dplyr::n_distinct(
        tf.gene.heatmap$gene
    )
)

message(
    "Relations TF-Gene après sélection du meilleur motif : ",
    nrow(tf.gene.cor)
)

message(
    "TF distincts : ",
    dplyr::n_distinct(
        tf.gene.cor$tf
    )
)

message(
    "Gènes distincts : ",
    dplyr::n_distinct(
        tf.gene.cor$gene
    )
)

message("Préparation de la matrice motif-peak...")

motif.matching <- objG@assays$peaks@motifs@data

if (
    ncol(motif.matching) !=
    ncol(motif_object)
) {

    stop(
        "Le nombre de motifs de motif.matching (",
        ncol(motif.matching),
        ") ne correspond pas au nombre de motifs de motif_object (",
        ncol(motif_object),
        ")."
    )
}

colnames(motif.matching) <-
    colnames(motif_object)

peaks.use <- intersect(
    unique(df.p2g$peak),
    rownames(motif.matching)
)


motif.matching <- motif.matching[
    peaks.use,
    ,
    drop = FALSE
]


message(
    "Peaks avec information motif : ",
    nrow(motif.matching)
)

message(
    "Motifs disponibles : ",
    ncol(motif.matching)
)

if (
    !all(
        grepl(
            "^MA[0-9]+\\.[0-9]+$",
            colnames(motif.matching)
        )
    )
) {

    stop(
        "Les colonnes de motif.matching ne sont pas ",
        "des identifiants JASPAR MAxxxx.x."
    )
}

selected_motif_map <- selected_tf_motifs_all %>%
    dplyr::select(
        motif_id,
        tf_label,
        tf
    ) %>%
    dplyr::distinct()


message(
    "Relations motif-TF retenues pour le GRN : ",
    nrow(selected_motif_map)
)


motif_ids_available <- intersect(
    selected_motif_map$motif_id,
    colnames(motif.matching)
)


selected_motif_map <- selected_motif_map %>%
    dplyr::filter(
        motif_id %in% motif_ids_available
    )


message(
    "Motifs sélectionnés présents dans motif.matching : ",
    dplyr::n_distinct(
        selected_motif_map$motif_id
    )
)

message(
    "Relations motif-TF conservées : ",
    nrow(selected_motif_map)
)


if (
    nrow(selected_motif_map) == 0
) {

    stop(
        "Aucun motif sélectionné n'est présent dans motif.matching."
    )
}


motif_columns <- match(
    selected_motif_map$motif_id,
    colnames(motif.matching)
)


if (
    anyNA(motif_columns)
) {

    stop(
        "Certains motif_id sélectionnés ne sont pas ",
        "présents dans motif.matching."
    )
}


motif.matching.selected <- motif.matching[
    ,
    motif_columns,
    drop = FALSE
]


colnames(
    motif.matching.selected
) <- selected_motif_map$tf


message(
    "Correspondance motif / TF pour le GRN : OK"
)

message(
    "Motifs uniques retenus : ",
    dplyr::n_distinct(
        selected_motif_map$motif_id
    )
)

message(
    "Relations motif-TF retenues : ",
    nrow(selected_motif_map)
)

message(
    "TF distincts retenus : ",
    dplyr::n_distinct(
        selected_motif_map$tf
    )
)


message("Construction du réseau de régulation...")

summ <- Matrix::summary(
    motif.matching.selected
)


df.p2m <- data.frame(

    peak =
        rownames(
            motif.matching.selected
        )[summ$i],

    tf =
        colnames(
            motif.matching.selected
        )[summ$j],

    is_bound =
        as.integer(
            summ$x
        ),

    stringsAsFactors = FALSE
)

if (
    nrow(df.p2m) == 0
) {

    stop(
        "Aucune relation Peak-TF n'a été trouvée."
    )
}


message(
    "Relations Peak-TF : ",
    nrow(df.p2m)
)

message(
    "TF distincts dans Peak-TF : ",
    dplyr::n_distinct(
        df.p2m$tf
    )
)

message(
    "Peaks distincts avec motif : ",
    dplyr::n_distinct(
        df.p2m$peak
    )
)

df.p2g_network <- df.p2g %>%
    dplyr::select(
        peak,
        gene
    ) %>%
    dplyr::distinct()

df.m2g <- dplyr::left_join(
    df.p2m,
    df.p2g_network,
    by = "peak"
) %>%
    dplyr::filter(
        !is.na(gene)
    ) %>%
    dplyr::group_by(
        tf,
        gene
    ) %>%
    dplyr::summarise(
        n_peaks = dplyr::n(),
        .groups = "drop"
    )


message(
    "Relations TF-Gene issues des motifs + Peak-to-Gene : ",
    nrow(df.m2g)
)

df.grn <- dplyr::left_join(
    df.m2g,
    tf.gene.cor,
    by = c(
        "tf",
        "gene"
    )
)


message(
    "Relations GRN avec corrélation : ",
    sum(
        is.finite(
            df.grn$correlation
        )
    )
)

df.grn_filt <- df.grn %>%
    dplyr::filter(
        is.finite(correlation),
        is.finite(p_value),
        is.finite(fdr)
    )


df.grn_filtered <- df.grn_filt %>%
    dplyr::filter(
        fdr < 0.0001
    )


message(
    "Relations GRN après FDR < 0.0001 : ",
    nrow(df.grn_filtered)
)


message(
    "TF distincts dans le GRN : ",
    dplyr::n_distinct(
        df.grn_filtered$tf
    )
)


message(
    "Gènes distincts dans le GRN : ",
    dplyr::n_distinct(
        df.grn_filtered$gene
    )
)

GRNHeatmap_adapted <- function(
    tf.gene.cor,
    tf.timepoint = NULL,
    gene.timepoint = NULL
){

    mat.cor <- tf.gene.cor %>%
        dplyr::select(
            tf,
            gene,
            correlation
        ) %>%
        dplyr::distinct() %>%
        tidyr::pivot_wider(
            names_from = tf,
            values_from = correlation
        ) %>%
        tibble::column_to_rownames(
            "gene"
        ) %>%
        as.matrix()

    if (!is.null(tf.timepoint)) {

        tf_present <- intersect(
            colnames(mat.cor),
            names(tf.timepoint)
        )

        mat.cor <- mat.cor[
            ,
            tf_present,
            drop = FALSE
        ]

        tf_order <- tf_present[
            order(
                tf.timepoint[tf_present],
                tf_present
            )
        ]

        mat.cor <- mat.cor[
            ,
            tf_order,
            drop = FALSE
        ]
    }

    if (!is.null(gene.timepoint)) {

        gene_present <- intersect(
            rownames(mat.cor),
            names(gene.timepoint)
        )

        mat.cor <- mat.cor[
            gene_present,
            ,
            drop = FALSE
        ]

        gene_order <- gene_present[
            order(
                gene.timepoint[gene_present],
                gene_present
            )
        ]

        mat.cor <- mat.cor[
            gene_order,
            ,
            drop = FALSE
        ]
    }

    column_ha <- NULL

    if (!is.null(tf.timepoint)) {

        tf_values <- tf.timepoint[
            colnames(mat.cor)
        ]

        names(tf_values) <-
            colnames(mat.cor)

        column_ha <-
            ComplexHeatmap::HeatmapAnnotation(
                TF_time = tf_values,
                annotation_name_side = "left"
            )
    }


    row_ha <- NULL

    if (!is.null(gene.timepoint)) {

        gene_values <- gene.timepoint[
            rownames(mat.cor)
        ]

        names(gene_values) <-
            rownames(mat.cor)

        row_ha <-
            ComplexHeatmap::rowAnnotation(
                gene_time = gene_values,
                show_legend = TRUE
            )
    }

    if (!is.null(tf.timepoint)) {

        if (
            !identical(
                names(tf_values),
                colnames(mat.cor)
            )
        ) {

            stop(
                "Annotation TF mal alignée avec la matrice."
            )
        }
    }


    if (!is.null(gene.timepoint)) {

        if (
            !identical(
                names(gene_values),
                rownames(mat.cor)
            )
        ) {

            stop(
                "Annotation gène mal alignée avec la matrice."
            )
        }
    }


    message(
        "Ordre temporel TF : OK"
    )

    message(
        "Ordre temporel gènes : OK"
    )

    ComplexHeatmap::Heatmap(

        mat.cor,

        name = "correlation",

        cluster_columns = FALSE,

        cluster_rows = FALSE,

        top_annotation =
            column_ha,

        left_annotation =
            row_ha,

        show_row_names = FALSE,

        show_column_names = TRUE,

        column_names_rot = 90,

        border = TRUE,

        row_title =
            "Genes ordered by pseudotime",

        column_title =
            "Transcription Factors",

        col =
            circlize::colorRamp2(
                c(-1, 0, 1),
                c(
                    "#2166AC",
                    "#F7F7F7",
                    "#B2182B"
                )
            ),

        use_raster = TRUE
    )
}

tf_time <- setNames(
    tf_time_points$time_point,
    tf_time_points$tf
)

gene_time <- setNames(
    df_gene_time_point$gene_time_point,
    df_gene_time_point$gene
)


tf_time <- tf_time[
    intersect(
        names(tf_time),
        unique(
            tf.gene.cor$tf
        )
    )
]


df.grn_heatmap_genes <- df.grn_filtered %>%

    dplyr::filter(
        is.finite(correlation),
        correlation > 0.25
    )


genes_heatmap <- unique(
    df.grn_heatmap_genes$gene
)


message(
    "Relations GRN finales utilisées pour la heatmap : ",
    nrow(
        df.grn_heatmap_genes
    )
)

message(
    "Gènes distincts dans le GRN final : ",
    length(
        genes_heatmap
    )
)

gene_time_heatmap <- gene_time[
    intersect(
        genes_heatmap,
        names(gene_time)
    )
]


genes_without_time <- setdiff(
    genes_heatmap,
    names(gene_time_heatmap)
)

message(
    "Gènes du GRN final avec gene_time : ",
    length(
        gene_time_heatmap
    )
)

message(
    "Gènes du GRN final sans gene_time : ",
    length(
        genes_without_time
    )
)


if (
    length(genes_without_time) > 0
) {

    warning(
        "Certains gènes du GRN final n'ont pas de gene_time."
    )

}

gene_time <- gene_time_heatmap


message(
    "Gènes retenus pour la heatmap : ",
    length(gene_time)
)

message(
    "TF retenus pour la heatmap : ",
    length(tf_time)
)

tf.gene.cor.heatmap <-
    tf.gene.cor %>%

    dplyr::filter(
        tf %in% names(tf_time),
        gene %in% names(gene_time),
        is.finite(correlation)
    ) %>%

    dplyr::select(
        tf,
        gene,
        correlation
    ) %>%

    dplyr::distinct()


expected_pairs <-
    length(tf_time) *
    length(gene_time)

observed_pairs <-
    nrow(
        tf.gene.cor.heatmap
    )

message(
    "Paires TF-Gene attendues pour la heatmap : ",
    expected_pairs
)

message(
    "Paires TF-Gene disponibles pour la heatmap : ",
    observed_pairs
)


if (
    observed_pairs != expected_pairs
) {

    warning(
        "La matrice TF-Gene de la heatmap n'est pas complète : ",
        observed_pairs,
        " / ",
        expected_pairs,
        " paires."
    )

} else {

    message(
        "Matrice TF-Gene complète : OK"
    )
}


ht_grn <- GRNHeatmap_adapted(

    tf.gene.cor =
        tf.gene.cor.heatmap,

    tf.timepoint =
        tf_time,

    gene.timepoint =
        gene_time
)

dir.create(
    "results/GRN",
    recursive = TRUE,
    showWarnings = FALSE
)

grn_heatmap_pdf <- file.path(
    "results/GRN",
    "GRN_heatmap_Cluster_G.pdf"
)

pdf(
    grn_heatmap_pdf,
    width = 12,
    height = 16,
    useDingbats = FALSE
)

ComplexHeatmap::draw(
    ht_grn,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
)

dev.off()

message(
    "Heatmap GRN sauvegardée : ",
    normalizePath(grn_heatmap_pdf)
)

ht_grn

message("")
message("========================================================")
message("F. ANALYSE DU RESEAU GRN")
message("========================================================")


gene_late_min <- 35
gene_late_max <- 73.5

tf_late_min <- 35
tf_late_max <- 73.5

gene_phase <- function(x) {

    dplyr::case_when(

        !is.finite(x) ~ NA_character_,

        x >= gene_late_min &
        x <= gene_late_max ~ "Late",

        TRUE ~ "Early"
    )
}


df_gene_phase <- data.frame(

    gene = names(gene_time),

    time_point =
        as.numeric(gene_time),

    time_phase =
        gene_phase(
            as.numeric(gene_time)
        ),

    stringsAsFactors = FALSE
)


message(
    "Gènes Early : ",
    sum(
        df_gene_phase$time_phase == "Early",
        na.rm = TRUE
    )
)

message(
    "Gènes Late : ",
    sum(
        df_gene_phase$time_phase == "Late",
        na.rm = TRUE
    )
)


tf_phase <- function(x) {

    dplyr::case_when(

        !is.finite(x) ~ NA_character_,

        x >= tf_late_min &
        x <= tf_late_max ~ "Late",

        TRUE ~ "Early"
    )
}


df_tf_phase <- data.frame(

    tf = names(tf_time),

    time_point =
        as.numeric(tf_time),

    time_phase =
        tf_phase(
            as.numeric(tf_time)
        ),

    stringsAsFactors = FALSE
)


message(
    "TF Early : ",
    sum(
        df_tf_phase$time_phase == "Early",
        na.rm = TRUE
    )
)

message(
    "TF Late : ",
    sum(
        df_tf_phase$time_phase == "Late",
        na.rm = TRUE
    )
)


message("")
message("TF Late :")

print(
    df_tf_phase %>%
        dplyr::filter(
            time_phase == "Late"
        ) %>%
        dplyr::arrange(
            time_point
        )
)

nodes_time <- dplyr::bind_rows(

    df_tf_phase %>%
        dplyr::transmute(
            name = tf,
            node_kind = "TF",
            time_point,
            time_phase
        ),

    df_gene_phase %>%
        dplyr::transmute(
            name = gene,
            node_kind = "gene",
            time_point,
            time_phase
        )

) %>%

    dplyr::distinct(
        name,
        .keep_all = TRUE
    )


message(
    "Nœuds temporels disponibles : ",
    nrow(nodes_time)
)


grn_fdr <- df.grn_filtered %>%

    dplyr::filter(

        is.finite(correlation),

        is.finite(fdr),

        fdr < 1e-4
    )


message(
    "Relations GRN après FDR < 1e-4 : ",
    nrow(grn_fdr)
)

grn_cor_cutoff <- 0.25


df.grn_final <- grn_fdr %>%

    dplyr::filter(

        correlation >
            grn_cor_cutoff
    )


message(
    "Relations GRN après FDR + correlation > ",
    grn_cor_cutoff,
    " : ",
    nrow(df.grn_final)
)


message(
    "TF distincts dans le GRN : ",
    length(
        unique(
            df.grn_final$tf
        )
    )
)

message(
    "Gènes distincts dans le GRN : ",
    length(
        unique(
            df.grn_final$gene
        )
    )
)

df.grn_final <- df.grn_final %>%

    dplyr::left_join(

        df_tf_phase %>%

            dplyr::rename(

                tf_time =
                    time_point,

                tf_phase =
                    time_phase
            ),

        by = "tf"
    ) %>%

    dplyr::left_join(

        df_gene_phase %>%

            dplyr::rename(

                gene_time =
                    time_point,

                gene_phase =
                    time_phase
            ),

        by = "gene"
    )


message(
    "Relations sans temps TF : ",
    sum(
        is.na(
            df.grn_final$tf_time
        )
    )
)

message(
    "Relations sans temps gène : ",
    sum(
        is.na(
            df.grn_final$gene_time
        )
    )
)

df.grn_final <- df.grn_final %>%

    dplyr::mutate(

        regulatory_phase =
            dplyr::case_when(

                tf_phase == "Early" &
                gene_phase == "Early" ~

                    "Early_TF_to_Early_gene",

                tf_phase == "Early" &
                gene_phase == "Late" ~

                    "Early_TF_to_Late_gene",

                tf_phase == "Late" &
                gene_phase == "Early" ~

                    "Late_TF_to_Early_gene",

                tf_phase == "Late" &
                gene_phase == "Late" ~

                    "Late_TF_to_Late_gene",

                TRUE ~
                    NA_character_
            )
    )


message("")
message("Types de relations temporelles :")

print(
    table(
        df.grn_final$regulatory_phase,
        useNA = "ifany"
    )
)


edges_all <- df.grn_final %>%

    dplyr::transmute(

        from =
            as.character(tf),

        to =
            as.character(gene),

        correlation =
            as.numeric(correlation),

        regulatory_score =
            as.numeric(correlation),

        t_stat =
            as.numeric(t_stat),

        p_value =
            as.numeric(p_value),

        fdr =
            as.numeric(fdr),

        regulatory_phase
    )


message(
    "Arêtes préparées : ",
    nrow(edges_all)
)

all_nodes <- unique(

    c(
        edges_all$from,
        edges_all$to
    )
)


nodes_all <- data.frame(

    name =
        all_nodes,

    stringsAsFactors = FALSE
)


g_grn <- igraph::graph_from_data_frame(

    d =
        edges_all,

    vertices =
        nodes_all,

    directed =
        TRUE
)


GRN_G <- tidygraph::as_tbl_graph(
    g_grn
)


message(
    "Nœuds du réseau : ",
    igraph::vcount(
        g_grn
    )
)

message(
    "Arêtes du réseau : ",
    igraph::ecount(
        g_grn
    )
)

node_table <- GRN_G %>%

    tidygraph::activate(
        nodes
    ) %>%

    tibble::as_tibble() %>%

    dplyr::left_join(

        nodes_time,

        by = "name"
    )


message(
    "Nœuds avec phase temporelle : ",
    sum(
        !is.na(
            node_table$time_phase
        )
    ),
    " / ",
    nrow(
        node_table
    )
)



GRN_G <- GRN_G %>%

    tidygraph::activate(
        nodes
    ) %>%

    dplyr::mutate(

        node_kind =
            node_table$node_kind,

        time_point =
            node_table$time_point,

        time_phase =
            node_table$time_phase
    )


GRN_G <- GRN_G %>%

    tidygraph::activate(
        nodes
    ) %>%

    dplyr::mutate(

        degree_in =
            tidygraph::centrality_degree(
                mode = "in",
                weights =
                    abs(
                        regulatory_score
                    )
            ),

        degree_out =
            tidygraph::centrality_degree(
                mode = "out",
                weights =
                    abs(
                        regulatory_score
                    )
            ),

        degree_total =
            tidygraph::centrality_degree(
                mode = "total"
            ),

        betweenness =
            tidygraph::centrality_betweenness(
                weights =
                    1 /
                    pmax(
                        abs(
                            regulatory_score
                        ),
                        1e-10
                    )
            ),

        closeness =
            tidygraph::centrality_closeness(),

        eigen =
            tidygraph::centrality_eigen(),

        pagerank =
            tidygraph::centrality_pagerank(
                weights =
                    abs(
                        regulatory_score
                    )
            )
    )


degree_in_cutoff <- quantile(

    GRN_G %>%

        tidygraph::activate(
            nodes
        ) %>%

        dplyr::pull(
            degree_in
        ),

    0.90,

    na.rm = TRUE
)


degree_out_cutoff <- quantile(

    GRN_G %>%

        tidygraph::activate(
            nodes
        ) %>%

        dplyr::pull(
            degree_out
        ),

    0.90,

    na.rm = TRUE
)


GRN_G <- GRN_G %>%

    tidygraph::activate(
        nodes
    ) %>%

    dplyr::mutate(

        hub_type =
            dplyr::case_when(

                degree_in >
                degree_in_cutoff &

                degree_out >
                degree_out_cutoff ~

                    "Hub integrateur",

                degree_out >
                degree_out_cutoff ~

                    "Hub regulateur",

                degree_in >
                degree_in_cutoff ~

                    "Hub cible",

                TRUE ~

                    "Standard"
            )
    )


message("")
message("Classification des hubs :")

print(

    GRN_G %>%

        tidygraph::activate(
            nodes
        ) %>%

        tibble::as_tibble() %>%

        dplyr::count(
            hub_type,
            time_phase,
            node_kind
        )
)

nodes_data_all <- GRN_G %>%

    tidygraph::activate(
        nodes
    ) %>%

    tibble::as_tibble()


library(ggplot2)
library(patchwork)
library(GGally)
library(scales)


p_dist <- nodes_data_all %>%

    dplyr::select(

        degree_in,
        degree_out,
        betweenness,
        pagerank
    ) %>%

    tidyr::pivot_longer(
        everything()
    ) %>%

    ggplot(
        aes(
            x = value
        )
    ) +

    geom_density() +

    facet_wrap(
        ~ name,
        scales = "free",
        ncol = 4
    ) +

    labs(

        title =
            "Distribution des métriques de centralité",

        x =
            "Valeur",

        y =
            "Densité"
    ) +

    theme_minimal()


p_dist

p_cor <- nodes_data_all %>%

    dplyr::select(

        degree_in,
        degree_out,
        betweenness,
        pagerank
    ) %>%

    GGally::ggpairs(

        upper =
            list(

                continuous =
                    GGally::wrap(
                        "cor",
                        size = 3
                    )
            ),

        lower =
            list(

                continuous =
                    GGally::wrap(
                        "points",
                        alpha = 0.3,
                        size = 0.5
                    )
            ),

        diag =
            list(

                continuous =
                    GGally::wrap(
                        "densityDiag",
                        alpha = 0.5
                    )
            )
    ) +

    theme_bw()


p_cor

top_hubs <- nodes_data_all %>%

    dplyr::filter(
        node_kind == "TF"
    ) %>%

    dplyr::arrange(
        dplyr::desc(
            degree_total
        )
    ) %>%

    dplyr::slice_head(
        n = 30
    )


p_hubs <- ggplot(

    top_hubs,

    aes(

        x =
            reorder(
                name,
                degree_total
            ),

        y =
            degree_total,

        fill =
            hub_type
    )
) +

    geom_col() +

    coord_flip() +

    labs(

        title =
            "Top 30 TF selon le degré total",

        x =
            "TF",

        y =
            "Degré total",

        fill =
            "Type de hub"
    ) +

    theme_minimal()


p_hubs

dir.create(

    "results/GRN/plots",

    recursive = TRUE,

    showWarnings = FALSE
)


ggsave(

    "results/GRN/plots/GRN_metric_distributions.pdf",

    p_dist,

    width = 12,

    height = 7
)


ggsave(

    "results/GRN/plots/GRN_metric_correlations.pdf",

    p_cor,

    width = 10,

    height = 10
)


ggsave(

    "results/GRN/plots/GRN_top_TF_hubs.pdf",

    p_hubs,

    width = 8,

    height = 10
)


nodes_gephi <- nodes_data_all %>%

    dplyr::transmute(

        Id =
            as.character(name),

        Label =
            as.character(name),

        node_kind =
            as.character(node_kind),

        time_phase =
            as.character(time_phase),

        time_point =
            as.numeric(time_point),

        degree_in =
            as.numeric(degree_in),

        degree_out =
            as.numeric(degree_out),

        degree_total =
            as.numeric(degree_total),

        betweenness =
            as.numeric(betweenness),

        pagerank =
            as.numeric(pagerank),

        hub_type =
            as.character(hub_type)
    )


node_names <- igraph::V(GRN_G)$name

edges_data_all <- GRN_G %>%

    tidygraph::activate(
        edges
    ) %>%

    tibble::as_tibble()

edges_gephi <- edges_data_all %>%

    dplyr::transmute(

        Source =
            node_names[as.integer(from)],

        Target =
            node_names[as.integer(to)],

        correlation =
            as.numeric(correlation),

        regulatory_score =
            as.numeric(correlation),

        t_stat =
            as.numeric(t_stat),

        p_value =
            as.numeric(p_value),

        fdr =
            as.numeric(fdr),

        regulatory_phase =
            as.character(regulatory_phase)
    )


node_ids <- unique(
    nodes_gephi$Id
)


missing_sources <- setdiff(
    unique(edges_gephi$Source),
    node_ids
)


missing_targets <- setdiff(
    unique(edges_gephi$Target),
    node_ids
)


message(
    "Nœuds Gephi : ",
    length(node_ids)
)

message(
    "Sources absentes des nodes : ",
    length(missing_sources)
)

message(
    "Targets absentes des nodes : ",
    length(missing_targets)
)


if (
    length(missing_sources) > 0
) {

    warning(
        "Sources absentes : ",
        paste(
            head(missing_sources, 20),
            collapse = ", "
        )
    )
}


if (
    length(missing_targets) > 0
) {

    warning(
        "Targets absentes : ",
        paste(
            head(missing_targets, 20),
            collapse = ", "
        )
    )
}

message(
    "Sources NA : ",
    sum(
        is.na(edges_gephi$Source)
    )
)

message(
    "Targets NA : ",
    sum(
        is.na(edges_gephi$Target)
    )
)


message(
    "Source exemple : ",
    edges_gephi$Source[1]
)

message(
    "Target exemple : ",
    edges_gephi$Target[1]
)

print(
    head(
        nodes_gephi
    )
)

print(
    head(
        edges_gephi
    )
)

readr::write_csv(

    nodes_gephi,

    "results/GRN/GRN_nodes_Gephi.csv"
)


readr::write_csv(

    edges_gephi,

    "results/GRN/GRN_edges_Gephi.csv"
)


message(
    "Nodes Gephi exportés : ",
    normalizePath(
        "results/GRN/GRN_nodes_Gephi.csv"
    )
)


message(
    "Edges Gephi exportées : ",
    normalizePath(
        "results/GRN/GRN_edges_Gephi.csv"
    )
)

message("")
message("========================================================")
message("DISTRIBUTION TEMPORELLE DES GENES")
message("========================================================")


gene_node_ids <- unique(
    as.integer(
        c(
            edges_data_all$from,
            edges_data_all$to
        )
    )
)


genes_grn_final <- nodes_data_all %>%

    dplyr::filter(
        row_number() %in% gene_node_ids,
        node_kind == "gene"
    ) %>%

    dplyr::distinct(
        name,
        .keep_all = TRUE
    ) %>%

    dplyr::pull(
        name
    )


message(
    "Gènes distincts dans le GRN final : ",
    length(genes_grn_final)
)


gene_time_map <- setNames(

    as.numeric(
        df_gene_time_point$gene_time_point
    ),

    as.character(
        df_gene_time_point$gene
    )

)


gene_time_grn <- gene_time_map[
    genes_grn_final
]

n_gene_with_time <- sum(
    is.finite(
        gene_time_grn
    )
)


n_gene_without_time <- sum(
    !is.finite(
        gene_time_grn
    )
)


message(
    "Gènes du GRN final avec gene_time : ",
    n_gene_with_time
)


message(
    "Gènes du GRN final sans gene_time : ",
    n_gene_without_time
)


if (
    n_gene_without_time > 0
) {

    message("")
    message(
        "Gènes sans gene_time :"
    )

    print(

        data.frame(
            gene =
                names(
                    gene_time_grn
                )[
                    !is.finite(
                        gene_time_grn
                    )
                ]
        )

    )

}


if (
    n_gene_with_time !=
    length(
        genes_grn_final
    )
) {

    warning(
        "Tous les gènes du GRN final ne possèdent pas de gene_time."
    )

}

gene_time_grn <- gene_time_grn[
    is.finite(
        gene_time_grn
    )
]


gene_time_grn <- gene_time_grn[
    order(
        gene_time_grn
    )
]



gene_phase_grn <- dplyr::case_when(

    gene_time_grn < 35 |
        gene_time_grn > 73.5
        ~ "Early",

    gene_time_grn >= 35 &
        gene_time_grn <= 73.5
        ~ "Late",

    TRUE
        ~ NA_character_

)


gene_time_distribution <- tibble::tibble(

    gene =
        names(
            gene_time_grn
        ),

    gene_time =
        as.numeric(
            gene_time_grn
        ),

    time_phase =
        gene_phase_grn

)


message("")
message(
    "Distribution temporelle :"
)

print(

    gene_time_distribution %>%

        dplyr::count(
            time_phase
        )

)


message("")
message(
    "Early : ",
    sum(
        gene_phase_grn == "Early",
        na.rm = TRUE
    )
)


message(
    "Late : ",
    sum(
        gene_phase_grn == "Late",
        na.rm = TRUE
    )
)


message("")
message(
    "Time rank min : ",
    min(
        gene_time_grn,
        na.rm = TRUE
    )
)


message(
    "Time rank max : ",
    max(
        gene_time_grn,
        na.rm = TRUE
    )
)


message("")
message(
    "Résumé du gene_time :"
)

print(
    summary(
        gene_time_grn
    )
)

p_gene_time_hist <- ggplot2::ggplot(

    gene_time_distribution,

    ggplot2::aes(
        x = gene_time
    )

) +

    ggplot2::geom_histogram(

        bins = 30,

        fill = "grey85",

        color = "grey35"

    ) +

    ggplot2::geom_vline(

        xintercept = 35,

        linetype = "dashed",

        linewidth = 0.8

    ) +

    ggplot2::geom_vline(

        xintercept = 73.5,

        linetype = "dashed",

        linewidth = 0.8

    ) +

    ggplot2::labs(

        title =
            "Distribution temporelle des 1647 gènes du GRN final",

        subtitle =
            "Early < 35 ou > 73.5 ; Late = 35–73.5",

        x =
            "Gene time rank",

        y =
            "Nombre de gènes"

    ) +

    ggplot2::theme_bw(

        base_size = 12

    )


p_gene_time_hist


p_gene_time_density <- ggplot2::ggplot(

    gene_time_distribution,

    ggplot2::aes(
        x = gene_time
    )

) +

    ggplot2::geom_density(

        linewidth = 0.9,

        fill = "grey85",

        alpha = 0.7

    ) +

    ggplot2::geom_vline(

        xintercept = 35,

        linetype = "dashed",

        linewidth = 0.8

    ) +

    ggplot2::geom_vline(

        xintercept = 73.5,

        linetype = "dashed",

        linewidth = 0.8

    ) +

    ggplot2::labs(

        title =
            "Densité des gènes du GRN final le long du time rank",

        subtitle =
            "1647 gènes",

        x =
            "Gene time rank",

        y =
            "Densité"

    ) +

    ggplot2::theme_bw(

        base_size = 12

    )


p_gene_time_density


p_gene_time_phase <- ggplot2::ggplot(

    gene_time_distribution,

    ggplot2::aes(

        x = gene_time,

        fill = time_phase

    )

) +

    ggplot2::geom_density(

        alpha = 0.45,

        linewidth = 0.8

    ) +

    ggplot2::geom_vline(

        xintercept = 35,

        linetype = "dashed",

        linewidth = 0.8

    ) +

    ggplot2::geom_vline(

        xintercept = 73.5,

        linetype = "dashed",

        linewidth = 0.8

    ) +

    ggplot2::labs(

        title =
            "Distribution temporelle des gènes du GRN final",

        subtitle =
            "Classification Early / Late selon le gene time",

        x =
            "Gene time rank",

        y =
            "Densité",

        fill =
            "Phase"

    ) +

    ggplot2::scale_fill_manual(

        values = c(

            Early = "#E76F51",

            Late = "#2A9D8F"

        )

    ) +

    ggplot2::theme_bw(

        base_size = 12

    ) +

    ggplot2::theme(

        legend.position = "top"

    )


p_gene_time_phase


dir.create(

    "results/GRN/plots/Genes",

    recursive = TRUE,

    showWarnings = FALSE

)


ggplot2::ggsave(

    "results/GRN/plots/Genes/GRN_gene_time_histogram.pdf",

    p_gene_time_hist,

    width = 8,

    height = 6

)


ggplot2::ggsave(

    "results/GRN/plots/Genes/GRN_gene_time_density.pdf",

    p_gene_time_density,

    width = 8,

    height = 6

)


ggplot2::ggsave(

    "results/GRN/plots/Genes/GRN_gene_time_density_phase.pdf",

    p_gene_time_phase,

    width = 8,

    height = 6

)


readr::write_csv(

    gene_time_distribution,

    "results/GRN/GRN_genes_time_distribution.csv"

)


message("")
message(
    "Plots de distribution temporelle exportés dans : ",
    normalizePath(
        "results/GRN/plots/Genes"
    )
)

message(
    "Table gene/time/phase exportée : ",
    normalizePath(
        "results/GRN/GRN_genes_time_distribution.csv"
    )
)

message("")
message("========================================================")
message("METRIQUES DES TF")
message("========================================================")

tf_plot_data <- nodes_data_all %>%
    dplyr::filter(
        node_kind == "TF"
    ) %>%
    dplyr::distinct(
        name,
        .keep_all = TRUE
    )


message(
    "TF dans le GRN final : ",
    nrow(tf_plot_data)
)

if (!"time_point" %in% names(tf_plot_data)) {

    stop(
        "La colonne 'time_point' est absente de nodes_data_all."
    )

}

if (!"time_phase" %in% names(tf_plot_data)) {

    stop(
        "La colonne 'time_phase' est absente de nodes_data_all."
    )

}


tf_without_time <- tf_plot_data %>%
    dplyr::filter(
        !is.finite(
            as.numeric(time_point)
        )
    )


if (nrow(tf_without_time) > 0) {

    stop(
        "TF du GRN final sans time_point : ",
        paste(
            tf_without_time$name,
            collapse = ", "
        )
    )

}


tf_without_phase <- tf_plot_data %>%
    dplyr::filter(
        is.na(time_phase)
    )


if (nrow(tf_without_phase) > 0) {

    stop(
        "TF du GRN final sans time_phase : ",
        paste(
            tf_without_phase$name,
            collapse = ", "
        )
    )

}


message(
    "TF avec phase temporelle : ",
    sum(
        !is.na(
            tf_plot_data$time_phase
        )
    ),
    " / ",
    nrow(tf_plot_data)
)

message("")
message("Distribution temporelle des TF :")

print(
    tf_plot_data %>%
        dplyr::count(
            time_phase
        )
)

message("")
message("TF du GRN final ordonnés par time_point :")

print(
    tf_plot_data %>%
        dplyr::select(
            name,
            time_point,
            time_phase,
            degree_in,
            degree_out,
            degree_total,
            betweenness,
            pagerank,
            dplyr::any_of("closeness"),
            dplyr::any_of("eigen"),
            hub_type
        ) %>%
        dplyr::arrange(
            time_point
        )
)

message("")
message("TF du GRN final ordonnés par time_point :")


print(

    tf_plot_data %>%

        dplyr::select(

            name,

            time_point,

            time_phase,

            degree_in,

            degree_out,

            degree_total,

            betweenness,

            pagerank,

            dplyr::any_of(
                "closeness"
            ),

            dplyr::any_of(
                "eigen"
            ),

            hub_type

        ) %>%

        dplyr::arrange(
            time_point
        )

)

message("")
message("========================================================")
message("CONTROLE DES METRIQUES TF")
message("========================================================")


tf_metrics <- c(

    "degree_total",

    "betweenness",

    "closeness",

    "eigen",

    "pagerank"

)


tf_metrics_available <- tf_metrics[

    tf_metrics %in%
        names(tf_plot_data)

]


tf_metrics_missing <- setdiff(

    tf_metrics,

    names(tf_plot_data)

)


message(
    "Métriques disponibles : ",
    paste(
        tf_metrics_available,
        collapse = ", "
    )
)


if (
    length(tf_metrics_missing) > 0
) {

    warning(
        "Métriques absentes : ",
        paste(
            tf_metrics_missing,
            collapse = ", "
        )
    )

}


plot_tf_metric <- function(

    metric,

    data = tf_plot_data,

    log_scale = FALSE

) {

    if (
        !metric %in% names(data)
    ) {

        stop(
            "Métrique absente de nodes_data_all : ",
            metric
        )

    }

    d <- data %>%

        dplyr::transmute(

            name,

            time_phase,

            value =
                as.numeric(
                    .data[[metric]]
                ),

            tie =
                dplyr::coalesce(
                    as.numeric(
                        degree_total
                    ),
                    0
                )

        ) %>%

        dplyr::filter(

            is.finite(
                value
            )

        )


    if (
        log_scale
    ) {

        d <- d %>%

            dplyr::filter(
                value > 0
            ) %>%

            dplyr::mutate(

                value_plot =
                    log10(
                        value
                    )

            )

    } else {

        d <- d %>%

            dplyr::mutate(

                value_plot =
                    value

            )

    }

    d <- d %>%

        dplyr::arrange(

            dplyr::desc(
                value
            ),

            dplyr::desc(
                tie
            ),

            name

        )


    d$name <- factor(

        d$name,

        levels =
            rev(
                d$name
            )

    )

    ggplot2::ggplot(

        d,

        ggplot2::aes(

            x =
                value_plot,

            y =
                name,

            color =
                time_phase

        )

    ) +

        ggplot2::geom_point(

            size = 3.5

        ) +

        ggplot2::labs(

            title =
                paste0(
                    metric,
                    " - TF du GRN final"
                ),

            x =
                if (
                    log_scale
                ) {

                    paste0(
                        metric,
                        " (log10)"
                    )

                } else {

                    metric

                },

            y =
                NULL,

            color =
                "Phase"

        ) +

        ggplot2::theme_minimal(

            base_size = 11

        ) +

        ggplot2::theme(

            legend.position =
                "top"

        ) +

        ggplot2::scale_color_manual(

            values = c(

                Early =
                    "#E76F51",

                Late =
                    "#2A9D8F",

                Both =
                    "#577590",

                Unknown =
                    "grey50"

            ),

            drop = FALSE

        )

}

p_tf_degree_total <-

    plot_tf_metric(

        metric =
            "degree_total",

        log_scale =
            TRUE

    )


p_tf_betweenness <-

    plot_tf_metric(

        metric =
            "betweenness",

        log_scale =
            TRUE

    )


p_tf_closeness <-

    plot_tf_metric(

        metric =
            "closeness",

        log_scale =
            TRUE

    )


p_tf_eigen <-

    plot_tf_metric(

        metric =
            "eigen",

        log_scale =
            FALSE

    )


p_tf_pagerank <-

    plot_tf_metric(

        metric =
            "pagerank",

        log_scale =
            TRUE

    )

p_tf_degree_total

p_tf_betweenness

p_tf_closeness

p_tf_eigen

p_tf_pagerank


dir.create(

    "results/GRN/plots/TF",

    recursive =
        TRUE,

    showWarnings =
        FALSE

)


ggplot2::ggsave(

    "results/GRN/plots/TF/TF_degree_total.pdf",

    p_tf_degree_total,

    width =
        8,

    height =
        10

)


ggplot2::ggsave(

    "results/GRN/plots/TF/TF_betweenness.pdf",

    p_tf_betweenness,

    width =
        8,

    height =
        10

)


ggplot2::ggsave(

    "results/GRN/plots/TF/TF_closeness.pdf",

    p_tf_closeness,

    width =
        8,

    height =
        10

)


ggplot2::ggsave(

    "results/GRN/plots/TF/TF_eigen.pdf",

    p_tf_eigen,

    width =
        8,

    height =
        10

)


ggplot2::ggsave(

    "results/GRN/plots/TF/TF_pagerank.pdf",

    p_tf_pagerank,

    width =
        8,

    height =
        10

)


message("")
message("Plots TF exportés dans : results/GRN/plots/TF/")


tf_grn_final <- nodes_data_all %>%

    dplyr::filter(
        node_kind == "TF"
    ) %>%

    dplyr::distinct(
        name,
        .keep_all = TRUE
    )


tf_final <- as.character(
    tf_grn_final$name
)


message(
    "TF dans le GRN final : ",
    length(tf_final)
)

tf_forced_late <- c(
    "Gli2",
    "Elf1",
    "Mafk",
    "Elk3",
    "Pitx1",
    "Elk1"
)


message("")
message("TF forcés Late présents dans le GRN :")

print(
    intersect(
        tf_forced_late,
        tf_final
    )
)


message("")
message("TF forcés Late absents du GRN :")

print(
    setdiff(
        tf_forced_late,
        tf_final
    )
)

tf_phase_final <- df_tf_phase %>%

    dplyr::transmute(

        tf =
            as.character(tf),

        time_point =
            as.numeric(time_point),

        time_phase_original =
            as.character(time_phase)

    ) %>%

    dplyr::filter(

        tf %in% tf_final

    ) %>%

    dplyr::distinct(

        tf,
        .keep_all = TRUE

    ) %>%

    dplyr::mutate(

        time_phase =
            dplyr::if_else(

                tf %in% tf_forced_late,

                "Late",

                time_phase_original

            )

    )

tf_early <- tf_phase_final %>%

    dplyr::filter(
        time_phase == "Early"
    ) %>%

    dplyr::pull(
        tf
    ) %>%

    unique()


tf_late <- tf_phase_final %>%

    dplyr::filter(
        time_phase == "Late"
    ) %>%

    dplyr::pull(
        tf
    ) %>%

    unique()


message("")
message("TF Early : ", length(tf_early))
message("TF Late  : ", length(tf_late))


message("")
message("TF Early :")

print(tf_early)


message("")
message("TF Late :")

print(tf_late)


if (
    length(intersect(tf_early, tf_late)) > 0
) {

    stop(
        "Un TF est simultanément Early et Late."
    )

}


if (
    length(setdiff(tf_final, c(tf_early, tf_late))) > 0
) {

    warning(
        "Certains TF du GRN final n'ont pas de classification temporelle."
    )

    print(
        setdiff(
            tf_final,
            c(tf_early, tf_late)
        )
    )

}


message("")
message("Contrôle des TF forcés Late :")

print(

    tf_phase_final %>%

        dplyr::filter(
            tf %in% tf_forced_late
        ) %>%

        dplyr::select(
            tf,
            time_point,
            time_phase_original,
            time_phase
        )

)


g_igraph <- tidygraph::as.igraph(
    GRN_G
)


message("")
message("Réseau utilisé pour les perturbations :")

message(
    "Nœuds : ",
    igraph::vcount(
        g_igraph
    )
)

message(
    "Arêtes : ",
    igraph::ecount(
        g_igraph
    )
)


perturb_tf_edges <- function(

    g,
    tf_name

) {

    if (
        !tf_name %in%
        igraph::V(g)$name
    ) {

        stop(
            "TF absent du réseau : ",
            tf_name
        )

    }


    v <- which(

        igraph::V(g)$name ==
            tf_name

    )

    outgoing_edges <- igraph::incident(

        g,

        v = v,

        mode = "out"

    )


    n_removed <- length(
        outgoing_edges
    )

    g_perturbed <- igraph::delete_edges(

        g,

        outgoing_edges

    )


    n_nodes <- igraph::vcount(
        g_perturbed
    )

    n_edges <- igraph::ecount(
        g_perturbed
    )


    n_components <- tryCatch(

        igraph::components(
            g_perturbed
        )$no,

        error = function(e)
            NA_integer_

    )


    isolated_nodes <- tryCatch(

        sum(
            igraph::degree(
                g_perturbed,
                mode = "all"
            ) == 0
        ),

        error = function(e)
            NA_integer_

    )


    density <- tryCatch(

        igraph::edge_density(
            g_perturbed
        ),

        error = function(e)
            NA_real_

    )


    tibble::tibble(

        tf =
            tf_name,

        n_edges_before =
            igraph::ecount(g),

        n_edges_removed =
            n_removed,

        n_edges_after =
            n_edges,

        fraction_edges_removed =
            n_removed /
            igraph::ecount(g),

        n_nodes =
            n_nodes,

        n_components =
            n_components,

        isolated_nodes =
            isolated_nodes,

        edge_density =
            density

    )

}

message("")
message("========================================================")
message("TEST PERTURBATION")
message("========================================================")


test_tf <- tf_late[1]


test_perturbation <- perturb_tf_edges(

    g_igraph,

    test_tf

)


print(
    test_perturbation
)


message("")
message("========================================================")
message("PERTURBATION TF EARLY")
message("========================================================")


perturbation_early <- purrr::map_dfr(

    tf_early,

    function(tf) {

        message(
            "Perturbation : ",
            tf
        )

        perturb_tf_edges(

            g_igraph,

            tf

        )

    }

)


message("")
message(
    "TF Early perturbés : ",
    nrow(
        perturbation_early
    )
)


message("")
message("========================================================")
message("PERTURBATION TF LATE")
message("========================================================")


perturbation_late <- purrr::map_dfr(

    tf_late,

    function(tf) {

        message(
            "Perturbation : ",
            tf
        )

        perturb_tf_edges(

            g_igraph,

            tf

        )

    }

)


message("")
message(
    "TF Late perturbés : ",
    nrow(
        perturbation_late
    )
)


perturbation_early <- perturbation_early %>%

    dplyr::arrange(

        dplyr::desc(
            n_edges_removed
        )

    ) %>%

    dplyr::mutate(

        impact_score =
            n_edges_removed

    )


perturbation_late <- perturbation_late %>%

    dplyr::arrange(

        dplyr::desc(
            n_edges_removed
        )

    ) %>%

    dplyr::mutate(

        impact_score =
            n_edges_removed

    )


top_TF_early <- perturbation_early %>%

    dplyr::slice_head(
        n = 20
    )


print(
    top_TF_early
)


top_TF_late <- perturbation_late %>%

    dplyr::slice_head(
        n = 20
    )



print(
    top_TF_late
)


p_perturbation_early <-

    ggplot2::ggplot(

        top_TF_early,

        ggplot2::aes(

            x =
                reorder(
                    tf,
                    n_edges_removed
                ),

            y =
                n_edges_removed

        )

    ) +

    ggplot2::geom_col() +

    ggplot2::coord_flip() +

    ggplot2::labs(

        title =
            "Top TF Early - perturbation des edges",

        x =
            "TF",

        y =
            "Edges sortantes supprimées"

    ) +

    ggplot2::theme_minimal()

p_perturbation_late <-

    ggplot2::ggplot(

        top_TF_late,

        ggplot2::aes(

            x =
                reorder(
                    tf,
                    n_edges_removed
                ),

            y =
                n_edges_removed

        )

    ) +

    ggplot2::geom_col() +

    ggplot2::coord_flip() +

    ggplot2::labs(

        title =
            "Top TF Late - perturbation des edges",

        x =
            "TF",

        y =
            "Edges sortantes supprimées"

    ) +

    ggplot2::theme_minimal()

p_perturbation_early

p_perturbation_late

dir.create(

    "results/GRN/perturbation",

    recursive = TRUE,

    showWarnings = FALSE

)


dir.create(

    "results/GRN/perturbation/plots",

    recursive = TRUE,

    showWarnings = FALSE

)

readr::write_csv(

    tf_phase_final,

    "results/GRN/perturbation/TF_phase.csv"

)


readr::write_csv(

    perturbation_early,

    "results/GRN/perturbation/TF_perturbation_Early.csv"

)


readr::write_csv(

    perturbation_late,

    "results/GRN/perturbation/TF_perturbation_Late.csv"

)

readr::write_csv(

    top_TF_early,

    "results/GRN/perturbation/Top_TF_Early.csv"

)


readr::write_csv(

    top_TF_late,

    "results/GRN/perturbation/Top_TF_Late.csv"

)


ggplot2::ggsave(

    "results/GRN/perturbation/plots/Top_TF_Early_edges.pdf",

    p_perturbation_early,

    width = 8,

    height = 8

)


ggplot2::ggsave(

    "results/GRN/perturbation/plots/Top_TF_Late_edges.pdf",

    p_perturbation_late,

    width = 8,

    height = 8

)

message("")
message("========================================================")
message("CONTROLE FINAL PERTURBATION")
message("========================================================")


message(
    "TF Early : ",
    length(tf_early)
)


message(
    "TF Late : ",
    length(tf_late)
)


message(
    "Perturbations Early : ",
    nrow(perturbation_early)
)


message(
    "Perturbations Late : ",
    nrow(perturbation_late)
)


message("")
message(
    "Edges du réseau initial : ",
    igraph::ecount(g_igraph)
)


message("")
message("Top 10 TF Early :")

print(

    perturbation_early %>%

        dplyr::slice_head(
            n = 10
        )

)


message("")
message("Top 10 TF Late :")

print(

    perturbation_late %>%

        dplyr::slice_head(
            n = 10
        )

)


message("")
message("Résultats exportés dans :")

message(
    normalizePath(
        "results/GRN/perturbation"
    )
)