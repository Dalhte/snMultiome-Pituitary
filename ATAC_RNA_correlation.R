
suppressPackageStartupMessages({

  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(pbapply)

})


options(
  stringsAsFactors = FALSE
)


cluster_name <- "Cluster_G"

trajectory_col <- paste0(
  "pt_slingshot2d_",
  cluster_name
)

trajectory_prefix <- "pt_slingshot2d_"

assay_atac <- "peaks"

assay_rna <- "SoupXRNA"



nbins <- 260

cells_per_bin_expected <- 10



smooth_k <- 9

spline_spar <- 0.45



dar_q_emp_threshold <- 0.05

deg_q_emp_threshold <- 0.05



promoter_upstream <- 500

promoter_downstream <- 500


obj_file <-
  "results/objects/PGintegrated73.reclustered.linkpeaks.trajectory.rds"


atac_perm_file <-
  "results/pseudotime_perm/Peaks/Cluster_G_Peak_perm.csv"


rna_perm_file <-
  "results/pseudotime_perm/RNA/Cluster_G_RNA_perm.csv"


tad_file <-
  "data/rat_boundaries_filtered.bed"


outdir <-
  "results/pseudotime_DE/ATAC_RNA_correlation"


dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


required_files <- c(
  obj_file,
  atac_perm_file,
  rna_perm_file,
  tad_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if(length(missing_files) > 0){

  stop(
    paste(
      "Fichiers manquants :",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )

}


cat(
  "\n========================================\n",
  "CHARGEMENT OBJET SEURAT\n",
  "========================================\n"
)

integrated73 <- readRDS(
  obj_file
)

cat(
  "Objet chargé.\n"
)


required_assays <- c(
  assay_atac,
  assay_rna
)

missing_assays <- setdiff(
  required_assays,
  names(integrated73)
)

if(length(missing_assays) > 0){

  stop(
    "Assays manquants : ",
    paste(
      missing_assays,
      collapse = ", "
    )
  )

}


meta <- integrated73[[]]


required_metadata <- c(
  "traj_cluster",
  trajectory_col
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(meta)
)

if(length(missing_metadata) > 0){

  stop(
    "Colonnes metadata manquantes : ",
    paste(
      missing_metadata,
      collapse = ", "
    )
  )

}


cells_keep <- rownames(meta)[

  meta$traj_cluster == cluster_name &

  is.finite(
    as.numeric(
      meta[[trajectory_col]]
    )
  )

]


if(length(cells_keep) == 0){

  stop(
    "Aucune cellule valide pour ",
    cluster_name
  )

}


cat(
  "Cluster : ",
  cluster_name,
  "\n",
  sep = ""
)

cat(
  "Cellules avec pseudotemps : ",
  length(cells_keep),
  "\n",
  sep = ""
)


obj <- subset(
  integrated73,
  cells = cells_keep
)



pt <- as.numeric(
  obj[[]][[trajectory_col]]
)


if(any(!is.finite(pt))){

  stop(
    "Pseudotemps non finis après sélection."
  )

}


pt <- (
  pt - min(pt)
) / (
  max(pt) - min(pt)
)



ord <- order(
  pt
)

bin <- integer(
  length(pt)
)

bin[ord] <-
  ceiling(
    seq_along(ord) *
      nbins /
      length(ord)
  )


cat(
  "Bins : ",
  nbins,
  "\n",
  sep = ""
)

cat(
  "Cellules/bin attendues : ~",
  round(
    length(pt) / nbins,
    2
  ),
  "\n",
  sep = ""
)



cat(
  "\n========================================\n",
  "LECTURE PERMUTATIONS\n",
  "========================================\n"
)


atac_perm <- fread(
  atac_perm_file
)


rna_perm <- fread(
  rna_perm_file
)



if(!all(
  c(
    "peak",
    "q_emp"
  ) %in%
  colnames(atac_perm)
)){

  stop(
    "Colonnes peak/q_emp absentes du fichier ATAC."
  )

}


if(!all(
  c(
    "gene",
    "q_emp"
  ) %in%
  colnames(rna_perm)
)){

  stop(
    "Colonnes gene/q_emp absentes du fichier RNA."
  )

}



DAR <- atac_perm %>%

  filter(
    !is.na(q_emp),
    q_emp < dar_q_emp_threshold
  ) %>%

  pull(
    peak
  ) %>%

  unique()


DEG <- rna_perm %>%

  filter(
    !is.na(q_emp),
    q_emp < deg_q_emp_threshold
  ) %>%

  pull(
    gene
  ) %>%

  unique()


cat(
  "DAR q_emp < ",
  dar_q_emp_threshold,
  " : ",
  length(DAR),
  "\n",
  sep = ""
)

cat(
  "DEG q_emp < ",
  deg_q_emp_threshold,
  " : ",
  length(DEG),
  "\n",
  sep = ""
)



DAR_perm <- unique(
  as.character(DAR)
)

DEG_perm <- unique(
  as.character(DEG)
)



peak_ids <- rownames(
  obj[[assay_atac]]
)

peak_ids_lower <- tolower(
  peak_ids
)

DAR_object <- peak_ids[
  peak_ids_lower %in% tolower(DAR_perm)
]


rna_features <- rownames(
  obj[[assay_rna]]
)

rna_features_lower <- tolower(
  rna_features
)

DEG_object <- rna_features[
  rna_features_lower %in% tolower(DEG_perm)
]

DEG_object_lower <- tolower(
  DEG_object
)


cat(
  "DAR présents dans objet : ",
  length(DAR_object),
  "\n",
  sep = ""
)

cat(
  "DEG présents dans objet : ",
  length(DEG_object),
  "\n",
  sep = ""
)


peak_gr_all <- granges(
  obj[[assay_atac]]
)

peak_gr_names <- rownames(
  obj[[assay_atac]]
)


DAR_idx <- which(
  peak_gr_names %in% DAR_object
)

peak_gr <- peak_gr_all[
  DAR_idx
]

peak_gr$peak_id <- peak_gr_names[
  DAR_idx
]



cat(
  "\n========================================\n",
  "CONSTRUCTION TADS\n",
  "========================================\n"
)

bnd_df <- fread(
  tad_file,
  col.names = c(
    "chr",
    "start",
    "end"
  )
)

bnd_df <- bnd_df %>%
  mutate(
    chr = as.character(chr),
    start = as.integer(start),
    end = as.integer(end)
  ) %>%
  arrange(
    chr,
    start
  )

build_tad_df <- function(df_chr){

  if(nrow(df_chr) < 2){

    return(
      data.frame(
        chr = character(),
        tad_start = integer(),
        tad_end = integer()
      )
    )

  }

  data.frame(

    chr =
      df_chr$chr[
        -nrow(df_chr)
      ],

    tad_start =
      df_chr$end[
        -nrow(df_chr)
      ] + 1L,

    tad_end =
      df_chr$start[
        -1
      ] - 1L,

    stringsAsFactors = FALSE

  )

}


tads_df <- bnd_df %>%

  group_split(
    chr,
    .keep = TRUE
  ) %>%

  lapply(
    build_tad_df
  ) %>%

  bind_rows() %>%

  mutate(
    TAD =
      paste0(
        chr,
        "_TAD_",
        row_number()
      )
  )


tad_gr <- makeGRangesFromDataFrame(

  tads_df,

  seqnames.field =
    "chr",

  start.field =
    "tad_start",

  end.field =
    "tad_end",

  keep.extra.columns =
    TRUE

)


cat(
  "\nCONTROLE TAD GRANGES\n"
)

cat(
  "TAD chromosomes :\n"
)

print(
  sort(
    unique(
      as.character(
        seqnames(tad_gr)
      )
    )
  )
)

cat(
  "Peak chromosomes :\n"
)

print(
  sort(
    unique(
      as.character(
        seqnames(peak_gr)
      )
    )
  )
)

cat(
  "TAD seqlevels :\n"
)

print(
  sort(
    seqlevels(tad_gr)
  )
)

cat(
  "Peak seqlevels :\n"
)

print(
  sort(
    seqlevels(peak_gr)
  )
)


cat(
  "\n========================================\n",
  "CONSTRUCTION TSS SIGNAC\n",
  "========================================\n"
)


ann <- Signac::Annotation(
  obj[["peaks"]]
)

ann_df <- as.data.frame(ann)

tss_df <- ann_df %>%
  mutate(
    strand = as.character(strand)
  ) %>%
  group_by(gene_name) %>%
  summarise(
    TSS_chr = first(as.character(seqnames)),
    strand = first(strand),
    TSS = if(first(strand) == "+"){
      min(start)
    } else {
      max(end)
    },
    .groups = "drop"
  )

tss_df$gene_original <- tss_df$gene_name

tss_df$gene <- tolower(
  tss_df$gene_name
)


rna_features <- rownames(
  obj[[assay_rna]]
)

rna_features_lower <- tolower(
  rna_features
)


tss_df <- tss_df %>%

  filter(
    gene %in% rna_features_lower
  )

gene_match <- match(
  tss_df$gene,
  rna_features_lower
)

tss_df$gene <- rna_features[
  gene_match
]

tss_df$gene_lower <- tolower(
  tss_df$gene
)

cat(
  "Genes RNA disponibles : ",
  length(rna_features),
  "\n",
  sep = ""
)

cat(
  "Genes RNA avec TSS : ",
  nrow(tss_df),
  "\n",
  sep = ""
)


tss_gr <- GRanges(

  seqnames =
    tss_df$TSS_chr,

  ranges =
    IRanges(

      start =
        tss_df$TSS,

      end =
        tss_df$TSS

    ),

  gene =
  tss_df$gene_original,

  strand =
    tss_df$strand

)


promoter_gr <- promoters(

  tss_gr,

  upstream =
    promoter_upstream,

  downstream =
    promoter_downstream

)

cat(
  "Peaks DAR avec coordonnées : ",
  length(peak_gr),
  "\n"
)

cat(
  "TADs : ",
  length(tad_gr),
  "\n"
)


ol_peaks_tads <- findOverlaps(

  peak_gr,

  tad_gr,

  ignore.strand = TRUE

)


peak_tad_df <- data.frame(

  peak =
    mcols(
      peak_gr
    )$peak_id[
      queryHits(
        ol_peaks_tads
      )
    ],

  TAD =
    mcols(
      tad_gr
    )$TAD[
      subjectHits(
        ol_peaks_tads
      )
    ]

) %>%

  distinct()


ol_genes_tad <- findOverlaps(

  tss_gr,

  tad_gr,

  ignore.strand = TRUE

)


gene_tad_df <- data.frame(

  gene =
    mcols(
      tss_gr
    )$gene[
      queryHits(
        ol_genes_tad
      )
    ],

  TAD =
    mcols(
      tad_gr
    )$TAD[
      subjectHits(
        ol_genes_tad
      )
    ]

) %>%

  distinct()



bnd_gr <- makeGRangesFromDataFrame(
  bnd_df
)

bnd_gr$idx <-
  seq_along(
    bnd_gr
  )


boundary_map <- tads_df %>%

  mutate(

    left_idx =
      as.integer(
        row_number()
      ),

    right_idx =
      as.integer(
        row_number() + 1
      )

  ) %>%

  pivot_longer(

    c(
      left_idx,
      right_idx
    ),

    values_to =
      "idx",

    names_to = NULL

  ) %>%

  filter(
    !is.na(idx)
  ) %>%

  dplyr::select(
    TAD,
    idx
  )


ol_genes_bnd <- findOverlaps(

  tss_gr,

  bnd_gr,

  ignore.strand = TRUE

)


if(length(ol_genes_bnd) > 0){

  genes_on_bnd <- data.frame(

    gene =
      mcols(
        tss_gr
      )$gene[
        queryHits(
          ol_genes_bnd
        )
      ],

    idx =
      mcols(
        bnd_gr
      )$idx[
        subjectHits(
          ol_genes_bnd
        )
      ]

  ) %>%

    left_join(
      boundary_map,
      by = "idx",
      relationship = "many-to-many"
    ) %>%

    dplyr::select(
      gene,
      TAD
    )


  gene_tad_df <-
    bind_rows(
      gene_tad_df,
      genes_on_bnd
    ) %>%

    distinct()

}


ol_prom <- findOverlaps(

  peak_gr,

  promoter_gr,

  ignore.strand = TRUE

)


if(length(ol_prom) > 0){

  prom_df <- data.frame(

    gene =
      mcols(
        promoter_gr
      )$gene[
        subjectHits(
          ol_prom
        )
      ],

    peak =
      mcols(
        peak_gr
      )$peak_id[
        queryHits(
          ol_prom
        )
      ],

peak_center =
  start(peak_gr)[
    queryHits(ol_prom)
  ] +
  (
    width(peak_gr)[
      queryHits(ol_prom)
    ] - 1
  ) / 2,

    tss_pos =
      start(
        promoter_gr
      )[
        subjectHits(
          ol_prom
        )
      ] +
      promoter_upstream

  )


  prom_df$dist <-
    abs(
      prom_df$peak_center -
      prom_df$tss_pos
    )


  promoter_links_df <-
    prom_df %>%

    group_by(
      gene
    ) %>%

    slice_min(
      dist,
      with_ties = FALSE
    ) %>%

    ungroup() %>%

    mutate(
      is_promoter = TRUE
    ) %>%

    dplyr::select(
      gene,
      peak,
      is_promoter
    )

} else {

  promoter_links_df <-
    data.frame(

      gene = character(),

      peak = character(),

      is_promoter = logical()

    )

}


raw_links <- peak_tad_df %>%

  inner_join(

    gene_tad_df,

    by = "TAD",

    relationship = "many-to-many"

  ) %>%

  dplyr::select(
    gene,
    peak,
    TAD
  )


raw_links <- raw_links %>%
  left_join(
    promoter_links_df,
    by = c("gene", "peak")
  ) %>%
  mutate(
    is_promoter = ifelse(
      is.na(is_promoter),
      FALSE,
      is_promoter
    )
  ) %>%
  distinct(
    gene,
    peak,
    .keep_all = TRUE
  )


cat(
  "\n========================================\n",
  "LIENS GENOMIQUES\n",
  "========================================\n"
)

cat(
  "Liens DAR-GENES : ",
  nrow(raw_links),
  "\n",
  sep = ""
)

cat(
  "Liens promoteurs : ",
  sum(raw_links$is_promoter),
  "\n",
  sep = ""
)



cat(
  "\n========================================\n",
  "LECTURE MATRICES OBSERVEES\n",
  "========================================\n"
)


DefaultAssay(
  obj
) <- assay_atac


atac_data <- GetAssayData(

  obj,

  assay =
    assay_atac,

  layer =
    "data"

)


DefaultAssay(
  obj
) <- assay_rna


rna_data <- GetAssayData(

  obj,

  assay =
    assay_rna,

  layer =
    "data"

)


if(!identical(
  colnames(atac_data),
  colnames(rna_data)
)){

  common_cells <- intersect(
    colnames(atac_data),
    colnames(rna_data)
  )

  if(length(common_cells) != ncol(obj)){

    stop(
      "Les cellules RNA et ATAC ne correspondent pas."
    )

  }

  atac_data <- atac_data[
    ,
    common_cells,
    drop = FALSE
  ]

  rna_data <- rna_data[
    ,
    common_cells,
    drop = FALSE
  ]

  pt <- pt[
    match(
      common_cells,
      colnames(obj)
    )
  ]

}


bin_order <- split(
  seq_along(pt),
  bin
)


bin_profile <- function(
    values,
    bin_index
){

  vapply(

    bin_index,

    function(idx){

      mean(
        values[idx],
        na.rm = TRUE
      )

    },

    numeric(1)

  )

}

adaptive_rollmean <- function(

    x,

    k = 9

){

  if(
    k %% 2 == 0
  ){

    stop(
      "smooth_k doit être impair."
    )

  }


  n <- length(x)

  half <- floor(
    k / 2
  )

  out <- numeric(
    n
  )


  for(i in seq_len(n)){

    left <- max(
      1,
      i - half
    )

    right <- min(
      n,
      i + half
    )

    values <- x[
      left:right
    ]

    out[i] <-
      mean(
        values,
        na.rm = TRUE
      )

  }


  out

}

process_profile <- function(

    values

){

  y_bin <- bin_profile(
    values,
    bin_order
  )


  y_roll <- adaptive_rollmean(
    y_bin,
    k = smooth_k
  )


  x <- seq(
    0,
    1,
    length.out = nbins
  )


  ok <- is.finite(
    y_roll
  )


  if(
    sum(ok) < 5
  ){

    return(
      NULL
    )

  }


  fit <- smooth.spline(

    x =
      x[ok],

    y =
      y_roll[ok],

    spar =
      spline_spar

  )


  pred <- predict(

    fit,

    x = x

  )


  as.numeric(
    pred$y
  )

}


cat(
  "\n========================================\n",
  "PROFILS ATAC\n",
  "========================================\n"
)


DAR_profiles <- lapply(

  DAR_object,

  function(pk){

    values <- as.numeric(
      atac_data[pk, ]
    )

    process_profile(
      values
    )

  }

)

names(
  DAR_profiles
) <- DAR_object

cat(
  "\n========================================\n",
  "PROFILS RNA\n",
  "========================================\n"
)


cat(
  "\n========================================\n",
  "PROFILS RNA — TOUS LES GENES\n",
  "========================================\n"
)

gene_profiles <- lapply(

  unique(tss_df$gene),

  function(g){

    values <- as.numeric(
      rna_data[g, ]
    )

    process_profile(
      values
    )

  }

)

names(
  gene_profiles
) <- unique(
  tss_df$gene
)

cat(
  "\n========================================\n",
  "CORRELATIONS DAR <-> GENES\n",
  "========================================\n"
)


calc_correlation <- function(

    peak,

    gene

){

  x <- DAR_profiles[[peak]]

y <- gene_profiles[[gene]]

  if(
    is.null(x) ||
    is.null(y)
  ){

    return(
      c(
        Correlation = NA_real_,
        TStat = NA_real_,
        Pval = NA_real_
      )
    )

  }


  if(
    sd(x) == 0 ||
    sd(y) == 0
  ){

    return(
      c(
        Correlation = NA_real_,
        TStat = NA_real_,
        Pval = NA_real_
      )
    )

  }


  r <- cor(
    x,
    y,
    method = "pearson"
  )


  if(
    !is.finite(r) ||
    abs(r) >= 1
  ){

    p <- if(
      abs(r) == 1
    ){
      0
    } else {
      NA_real_
    }

  } else {

    t_stat <-
      r *
      sqrt(
        (
          length(x) - 2
        ) /
        (
          1 - r^2
        )
      )

    p <-
      2 *
      pt(
        -abs(t_stat),
        df =
          length(x) - 2
      )

  }


  if(
    abs(r) >= 1
  ){

    t_stat <- NA_real_

  }


  c(

    Correlation =
      r,

    TStat =
      t_stat,

    Pval =
      p

  )

}


res <- pbmapply(

  calc_correlation,

  raw_links$peak,

  raw_links$gene

)


res <- t(
  res
)


colnames(res) <- c(
  "Correlation",
  "TStat",
  "Pval"
)


raw_links <- raw_links %>%

  mutate(

    Correlation =
      as.numeric(
        res[
          ,
          "Correlation"
        ]
      ),

    TStat =
      as.numeric(
        res[
          ,
          "TStat"
        ]
      ),

    Pval =
      as.numeric(
        res[
          ,
          "Pval"
        ]
      )

  )


raw_links$FDR <-
  p.adjust(
    raw_links$Pval,
    method = "BH"
  )


peak_centers <- start(
  peak_gr_all
) +
  (
    width(
      peak_gr_all
    ) - 1
  ) / 2


names(
  peak_centers
) <-
  rownames(
    obj[[assay_atac]]
  )


tss_pos <- tss_df$TSS

names(tss_pos) <-
  tss_df$gene_original

raw_links$distance <-
  abs(

    peak_centers[
      raw_links$peak
    ] -

    tss_pos[
      raw_links$gene
    ]

  )


raw_file <- file.path(

  outdir,

  paste0(
    cluster_name,
    "_Peak_Gene_TAD_links_raw_correlation_spline.csv"
  )

)


fwrite(

  raw_links,

  raw_file

)

links_activator <- raw_links %>%

filter(
  !is.na(Correlation),
  !is.na(FDR),
  FDR < 0.001,
  Correlation > 0
) %>%

  dplyr::select(

    gene,
    peak,
    TAD,
    is_promoter,
    Correlation,
    TStat,
    Pval,
    FDR,
    distance

  )


links_silencer <- raw_links %>%

filter(
  !is.na(Correlation),
  !is.na(FDR),
  FDR < 0.001,
  Correlation < 0
) %>%

  dplyr::select(

    gene,
    peak,
    TAD,
    is_promoter,
    Correlation,
    TStat,
    Pval,
    FDR,
    distance

  )


activator_file <- file.path(

  outdir,

  paste0(
    cluster_name,
    "_Peak_Gene_TAD_links_activator.csv"
  )

)


silencer_file <- file.path(

  outdir,

  paste0(
    cluster_name,
    "_Peak_Gene_TAD_links_silencer.csv"
  )

)


fwrite(
  links_activator,
  activator_file
)


fwrite(
  links_silencer,
  silencer_file
)


cat(
  "\n========================================\n",
  "15 — ATAC / RNA CORRELATION TERMINEE\n",
  "========================================\n"
)

cat(
  "Cluster : ",
  cluster_name,
  "\n",
  sep = ""
)

cat(
  "Cellules : ",
  length(cells_keep),
  "\n",
  sep = ""
)

cat(
  "Bins : ",
  nbins,
  "\n",
  sep = ""
)

cat(
  "DAR : ",
  length(DAR),
  "\n",
  sep = ""
)

cat(
  "DEG : ",
  length(DEG),
  "\n",
  sep = ""
)

cat(
  "Liens DAR-GENES : ",
  nrow(raw_links),
  "\n",
  sep = ""
)

cat(
  "Liens positifs : ",
  nrow(links_activator),
  "\n",
  sep = ""
)

cat(
  "Liens négatifs : ",
  nrow(links_silencer),
  "\n",
  sep = ""
)

cat(
  "Smoothing : adaptive rollmean k = ",
  smooth_k,
  "\n",
  sep = ""
)

cat(
  "Spline spar = ",
  spline_spar,
  "\n",
  sep = ""
)

cat(
  "DAR q_emp < ",
  dar_q_emp_threshold,
  "\n",
  sep = ""
)

cat(
  "DEG q_emp < ",
  deg_q_emp_threshold,
  "\n",
  sep = ""
)

cat(
  "\nRAW : ",
  raw_file,
  "\n",
  sep = ""
)

cat(
  "ACTIVATOR : ",
  activator_file,
  "\n",
  sep = ""
)

cat(
  "SILENCER : ",
  silencer_file,
  "\n",
  sep = ""
)

cat(
  "========================================\n"
)