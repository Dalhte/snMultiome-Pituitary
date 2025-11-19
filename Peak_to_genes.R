# ---------------------------------------------------------------------------
# peak_to_gene_pipeline.R — robust to Seurat v5 (only 'peaks' present)
# ---------------------------------------------------------------------------
# What it does
#   1) Build peaks/TADs/promoters
#   2) Assemble candidate peak–gene links (TAD co-membership + promoter)
#   3) Correlate ATAC vs expression (raw or spline-smoothed)
#   4) Estimate empirical r-thresholds by circular permutation
#   5) Export raw links + activators + silencers
#
# Inputs
#   - Seurat object in memory: integrated73 (preferred) or integrated
#     * DefaultAssay must be or include a ChromatinAssay named "peaks"
#   - Expression matrix (genes × cells):
#       (A) If an RNA assay exists → auto-use (SoupXRNA or RNA)
#       (B) Else: provide CSV path with genes in rows, cell barcodes as columns
#   - BED files: TAD boundaries (3 cols), unique TSS BED6
#
# Outputs (under out_dir)
#   - Peak_Gene_TAD_links_raw_correlation[_spline].csv
#   - Peak_Gene_TAD_links_activator_FDR001_empirical[_spline].csv
#   - Peak_Gene_TAD_links_silencer_FDR001_empirical[_spline].csv
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(pbapply)
  library(readr)
  library(stringr)
})

# ---------------------------
# 0) Configuration (edit)
# ---------------------------
Config <- list(
  # BED inputs
  tad_boundaries_bed = "C:/Users/Charles/Desktop/Charles/SnMultiome/Data/rat_boundaries_filtered.bed",
  tss_unique_bed     = "C:/Users/Charles/Desktop/Charles/SnMultiome/Test final scripts/gene_TSS_unique_allGenes.bed",
  
  # Optional expression CSV (if no RNA assay in object)
  #   - rows: gene symbols (case-insensitive)
  #   - columns: cell barcodes matching the ATAC matrix columns
  expr_csv           = NA_character_,  # e.g. "C:/path/expression_matrix.csv"
  
  # Output directory (use absolute or project-root relative)
  out_dir            = "C:/Users/Charles/Desktop/Charles/SnMultiome/results/peak_to_gene/",
  
  # Parameters
  promoter_window_bp = 500L,           # TSS ± window
  blacklist_prefixes = c("^loc", "^newgene", "^rgd"),
  use_spline         = TRUE,           # smooth then correlate
  spline_spar        = 0.45,           # smooth.spline spar
  fdr_threshold      = 0.01,           # BH FDR cutoff for initial filter
  n_perm             = 1000,           # circular permutations for r cutoffs
  max_perm_rows      = 50000,          # cap #non-promoter rows sampled
  seed               = 123             # reproducibility
)

# Create output dir (recursive)
dir.create(Config$out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(Config$out_dir))

# ---------------------------
# 1) Inputs (robust grab)
# ---------------------------
seurat_obj <- if (exists("integrated")) integrated else if (exists("integrated_rds")) integrated_rds else NULL
stopifnot(inherits(seurat_obj, "Seurat"))

# --- ATAC matrix (from ChromatinAssay 'peaks') ---
if (!is.null(seurat_obj[["peaks"]])) {
  # prefer 'data' slot; fall back to 'counts'
  slot_use <- if ("data" %in% slotNames(seurat_obj[["peaks"]])) "data" else "counts"
  atac_mat <- GetAssayData(seurat_obj, assay = "peaks", slot = slot_use)
} else {
  stop("ChromatinAssay 'peaks' not found in the object.")
}

# --- Expression matrix ---
rna_assay <- if ("SoupXRNA" %in% names(seurat_obj@assays)) "SoupXRNA" else if ("RNA" %in% names(seurat_obj@assays)) "RNA" else NA_character_

if (!is.na(rna_assay)) {
  expr_mat <- GetAssayData(seurat_obj, assay = rna_assay, slot = "data")
} else {
  # No RNA assay → require external CSV
  if (is.na(Config$expr_csv) || !file.exists(Config$expr_csv)) {
    stop("No RNA assay in object and Config$expr_csv is missing or not found. Set Config$expr_csv to a valid gene×cell CSV.")
  }
  expr_df <- suppressMessages(readr::read_csv(Config$expr_csv, show_col_types = FALSE))
  # First column must be gene names if readr collapsed rownames; try to detect
  if (!any(tolower(names(expr_df)) %in% c("gene","genes","symbol","rownames"))) {
    # assume first column are rownames
    rn <- expr_df[[1]]
    expr_df <- expr_df[,-1, drop = FALSE]
    rownames(expr_df) <- rn
  } else {
    gene_col <- names(expr_df)[match(TRUE, tolower(names(expr_df)) %in% c("gene","genes","symbol","rownames"))]
    rn <- expr_df[[gene_col]]
    expr_df[[gene_col]] <- NULL
    rownames(expr_df) <- rn
  }
  expr_mat <- as.matrix(expr_df)
  rm(expr_df)
}

# --- Harmonize columns (cells) ---
common_cells <- intersect(colnames(atac_mat), colnames(expr_mat))
if (length(common_cells) < 10)
  stop("Too few shared cells between ATAC and expression. Check column names/barcodes.")
atac_mat <- atac_mat[, common_cells, drop = FALSE]
expr_mat <- expr_mat[, common_cells, drop = FALSE]

# --- Clean gene/peak identifiers (lowercase for joins) ---
rownames(expr_mat) <- tolower(rownames(expr_mat))
peak_ids <- tolower(rownames(atac_mat))

# ---------------------------
# 2) Genomic objects
# ---------------------------
# Peaks → GRanges (expects peak_id format 'chr-start-end' or 'chr:start-end')
parse_peak <- function(x) {
  x <- gsub(":", "-", x)
  data.frame(
    chr   = sub("-.*", "", x),
    start = as.integer(sub(".*?-([0-9]+)-.*", "\\1", x)),
    end   = as.integer(sub(".*-([0-9]+)$", "\\1", x)),
    stringsAsFactors = FALSE
  )
}
pk_df <- parse_peak(peak_ids)
if (any(!is.finite(pk_df$start) | !is.finite(pk_df$end)))
  stop("Peak IDs must be 'chr-start-end' (or 'chr:start-end'). Found non-parsable IDs.")
peak_gr <- GRanges(seqnames = pk_df$chr,
                   ranges   = IRanges(start = pk_df$start, end = pk_df$end),
                   peak_id  = peak_ids)

# TADs from boundaries (3 cols: chr start end)
message("Reading TAD boundaries: ", Config$tad_boundaries_bed)
bnd_df <- fread(Config$tad_boundaries_bed, col.names = c("chr","start","end")) |>
  mutate(chr = tolower(chr)) |>
  arrange(chr, as.integer(start), as.integer(end)) |>
  mutate(bnd_id = row_number())

build_tads <- function(df_chr) {
  if (nrow(df_chr) < 2) return(NULL)
  tibble(
    chr          = df_chr$chr[-nrow(df_chr)],
    tad_start    = as.integer(df_chr$end[-nrow(df_chr)]) + 1L,
    tad_end      = as.integer(df_chr$start[-1]) - 1L,
    left_bnd_id  = df_chr$bnd_id[-nrow(df_chr)],
    right_bnd_id = df_chr$bnd_id[-1]
  ) |>
    filter(tad_end >= tad_start) |>
    mutate(TAD = paste0(chr, ":TAD:", row_number()))
}
tads_df <- bnd_df |>
  group_by(chr) |>
  group_split(keep = TRUE) |>
  lapply(build_tads) |>
  bind_rows()
if (nrow(tads_df) == 0) stop("No TADs could be constructed from boundaries.")

tad_gr <- makeGRangesFromDataFrame(
  tads_df, seqnames.field = "chr",
  start.field = "tad_start", end.field = "tad_end",
  keep.extra.columns = TRUE
)

# Unique TSS BED6
message("Reading unique TSS BED: ", Config$tss_unique_bed)
tss_bed <- fread(
  Config$tss_unique_bed,
  col.names = c("chr","tss_start","tss_end","gene","score","strand")
) |>
  mutate(across(c(chr,gene,strand), tolower))

# Filter genes to those present in expr_mat and not blacklisted
bad_pat <- paste(Config$blacklist_prefixes, collapse = "|")
valid_genes <- rownames(expr_mat)[!grepl(bad_pat, rownames(expr_mat))]
tss_bed <- tss_bed |> filter(gene %in% valid_genes)

tss_gr <- GRanges(
  seqnames = tss_bed$chr,
  ranges   = IRanges(start = tss_bed$tss_start, end = tss_bed$tss_end),
  gene     = tss_bed$gene,
  strand   = "*"
)

promoter_gr <- promoters(
  tss_gr,
  upstream   = Config$promoter_window_bp,
  downstream = Config$promoter_window_bp
)

# ---------------------------
# 3) Assignments (peaks↔TAD; genes↔TAD; promoters)
# ---------------------------
# Peaks → TAD
ol_peaks_tads <- findOverlaps(peak_gr, tad_gr, ignore.strand = TRUE)
peak_tad_df <- tibble(
  peak = mcols(peak_gr)$peak_id[queryHits(ol_peaks_tads)],
  TAD  = mcols(tad_gr)$TAD[subjectHits(ol_peaks_tads)]
) |> distinct()

# Genes → TAD (by TSS overlap)
ol_genes_tad <- findOverlaps(tss_gr, tad_gr, ignore.strand = TRUE)
gene_tad_df <- tibble(
  gene = mcols(tss_gr)$gene[queryHits(ol_genes_tad)],
  TAD  = mcols(tad_gr)$TAD[subjectHits(ol_genes_tad)]
) |> distinct()

# Genes overlapping boundaries → link to adjacent TADs
bnd_gr <- makeGRangesFromDataFrame(bnd_df, keep.extra.columns = TRUE)
ol_genes_bnd <- findOverlaps(tss_gr, bnd_gr, ignore.strand = TRUE)
if (length(ol_genes_bnd) > 0) {
  genes_on_bnd <- tibble(
    gene   = mcols(tss_gr)$gene[queryHits(ol_genes_bnd)],
    bnd_id = mcols(bnd_gr)$bnd_id[subjectHits(ol_genes_bnd)]
  ) |>
    left_join(
      tads_df |>
        select(TAD, left_bnd_id, right_bnd_id) |>
        pivot_longer(c(left_bnd_id, right_bnd_id), values_to = "bnd_id") |>
        select(TAD, bnd_id) |>
        distinct(),
      by = "bnd_id"
    ) |>
    filter(!is.na(TAD)) |>
    select(gene, TAD) |>
    distinct()
  gene_tad_df <- bind_rows(gene_tad_df, genes_on_bnd) |> distinct()
}

# Promoter links: nearest peak per gene within promoter window
ol_prom <- findOverlaps(peak_gr, promoter_gr, ignore.strand = TRUE)
promoter_links_df <- if (length(ol_prom) > 0) {
  prom_df <- tibble(
    gene        = mcols(promoter_gr)$gene[subjectHits(ol_prom)],
    peak        = mcols(peak_gr)$peak_id[queryHits(ol_prom)],
    peak_center = start(peak_gr)[queryHits(ol_prom)] + (width(peak_gr)[queryHits(ol_prom)] - 1)/2,
    tss_pos     = start(promoter_gr)[subjectHits(ol_prom)] + Config$promoter_window_bp
  ) |>
    mutate(dist = abs(peak_center - tss_pos)) |>
    group_by(gene) |>
    slice_min(dist, with_ties = FALSE) |>
    ungroup() |>
    mutate(is_promoter = TRUE) |>
    select(gene, peak, is_promoter)
} else tibble(gene = character(), peak = character(), is_promoter = logical())

# ---------------------------------------------------------------------------
# 4) Technical columns (indices, distances)
# ---------------------------------------------------------------------------
# Build raw peak–gene–TAD table using the same logic as the original script:
#   1) Peaks and genes must share a TAD
#   2) If a peak is a promoter for at least one gene, remove that peak from all
#      TAD-based links and re-add only the promoter links.

raw_links <- peak_tad_df |>
  inner_join(gene_tad_df, by = "TAD", relationship = "many-to-many") |>
  dplyr::select(gene, peak, TAD)

# Remove all TAD links involving peaks that are used as promoters
raw_links <- raw_links |>
  dplyr::filter(!peak %in% promoter_links_df$peak)

# Add promoter links (they may not have a TAD assigned)
raw_links <- dplyr::bind_rows(raw_links, promoter_links_df) |>
  dplyr::mutate(
    is_promoter = ifelse(is.na(is_promoter), FALSE, is_promoter)
  )

# Compute gene/peak indices in the matrices and genomic distance (peak center ↔ TSS)
peak_centers <- start(peak_gr) + width(peak_gr) / 2
names(peak_centers) <- peak_ids

tss_pos <- start(tss_gr) + 0.5
names(tss_pos) <- mcols(tss_gr)$gene

raw_links <- raw_links |>
  mutate(
    gene_idx = match(gene, rownames(expr_mat)),
    peak_idx = match(peak, peak_ids),
    distance = abs(peak_centers[peak] - tss_pos[gene])
  ) |>
  # Keep only rows that map correctly to both matrices
  filter(!is.na(gene_idx) & !is.na(peak_idx))

# ---------------------------------------------------------------------------
# 5) Correlations (with optional spline smoothing)
# ---------------------------------------------------------------------------

apply_spline <- function(x, spar = Config$spline_spar) {
  t <- seq_along(x)
  y <- tryCatch({
    fit <- smooth.spline(t, x, spar = spar)
    as.numeric(predict(fit, t)$y)
  }, error = function(e) as.numeric(x))
  y
}

calc_cor <- function(pk_i, gn_i, use_spline = Config$use_spline) {
  x <- atac_mat[pk_i, ]
  y <- expr_mat[gn_i, ]
  if (all(is.na(x)) || all(is.na(y))) return(c(Correlation = NA, TStat = NA, Pval = NA))
  if (sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) return(c(Correlation = NA, TStat = NA, Pval = NA))
  if (use_spline) {
    x <- apply_spline(x)
    y <- apply_spline(y)
  }
  r <- suppressWarnings(cor(x, y, use = "complete.obs"))
  n <- sum(is.finite(x) & is.finite(y))
  if (!is.finite(r) || n < 3) return(c(Correlation = NA, TStat = NA, Pval = NA))
  t <- r * sqrt((n - 2) / max(1e-8, 1 - r^2))
  p <- 2 * pt(-abs(t), df = n - 2)
  c(Correlation = r, TStat = t, Pval = p)
}

message("Computing peak–gene correlations (", ifelse(Config$use_spline, "spline", "raw"), ") ...")
res_mat <- pbmapply(
  FUN = function(peak_idx, gene_idx) calc_cor(peak_idx, gene_idx),
  raw_links$peak_idx, raw_links$gene_idx
)
res_mat <- t(res_mat)
colnames(res_mat) <- c("Correlation", "TStat", "Pval")

raw_links <- bind_cols(raw_links, as.data.frame(res_mat)) |>
  mutate(FDR = p.adjust(Pval, method = "BH"))

raw_out <- file.path(
  Config$out_dir,
  sprintf(
    "Peak_Gene_TAD_links_raw_correlation%s.csv",
    ifelse(Config$use_spline, "_spline", "")
  )
)
write_csv(raw_links, raw_out)
message("✔ Raw links written: ", raw_out)

# ---------------------------------------------------------------------------
# 6) Empirical correlation thresholds + activator/silencer exports
# ---------------------------------------------------------------------------

set.seed(Config$seed)

# Normalize promoter flag and enforce lowercase identifiers
links_df <- raw_links |>
  mutate(
    is_promoter = !is.na(is_promoter) & is_promoter,
    peak = tolower(peak),
    gene = tolower(gene)
  )

fdr_thr <- Config$fdr_threshold

# Significant links:
#   - all promoter links are kept
#   - non-promoter links must pass the FDR threshold
links_sig_all <- links_df |>
  filter(is_promoter | (!is.na(FDR) & FDR < fdr_thr))

# Non-promoter links used to estimate the null distribution
non_prom <- links_sig_all |>
  filter(!is_promoter)

if (nrow(non_prom) == 0) {
  stop("No non-promoter significant links available to estimate empirical thresholds.")
}

# Sample a subset for permutations to limit runtime
sample_idx <- sample(
  seq_len(nrow(non_prom)),
  min(Config$max_perm_rows, nrow(non_prom))
)
perm_tbl <- non_prom[sample_idx, ]

# Permutation correlation:
#   - recompute correlations on *raw* signals (no spline),
#   - apply a random circular shift on expression profiles
#     to break the true association while preserving autocorrelation.
perm_cor_once <- function(i) {
  pk_i <- perm_tbl$peak_idx[i]
  gn_i <- perm_tbl$gene_idx[i]
  x <- atac_mat[pk_i, ]
  y <- expr_mat[gn_i, ]
  if (length(x) <= 1L) return(NA_real_)
  shift <- sample.int(length(x) - 1L, 1)
  y <- c(y[(shift + 1L):length(y)], y[1L:shift])
  suppressWarnings(cor(x, y, use = "complete.obs"))
}

nperm <- Config$n_perm
perm_vals <- replicate(
  nperm,
  perm_cor_once(sample.int(nrow(perm_tbl), 1)),
  simplify = TRUE
)

thr_pos <- as.numeric(quantile(perm_vals, 0.95, na.rm = TRUE))  # activators
thr_neg <- as.numeric(quantile(perm_vals, 0.05, na.rm = TRUE))  # silencers
message(sprintf("Empirical thresholds (raw r): r_pos = %.3f | r_neg = %.3f", thr_pos, thr_neg))

# Final activator and silencer sets:
#   - all promoter links are always kept as activators
#   - non-promoter activators: r > r_pos
#   - silencers: r < r_neg
links_activator <- links_sig_all |>
  filter(is_promoter | (!is_promoter & Correlation > thr_pos)) |>
  dplyr::select(
    gene, peak, TAD, is_promoter,
    Correlation, TStat, Pval, FDR, distance
  )

links_silencer <- links_sig_all |>
  filter(!is_promoter & Correlation < thr_neg) |>
  dplyr::select(
    gene, peak, TAD, is_promoter,
    Correlation, TStat, Pval, FDR, distance
  )

act_out <- file.path(
  Config$out_dir,
  sprintf(
    "Peak_Gene_TAD_links_activator_FDR001_empirical%s.csv",
    ifelse(Config$use_spline, "_spline", "")
  )
)
sil_out <- file.path(
  Config$out_dir,
  sprintf(
    "Peak_Gene_TAD_links_silencer_FDR001_empirical%s.csv",
    ifelse(Config$use_spline, "_spline", "")
  )
)

fwrite(links_activator, act_out)
fwrite(links_silencer,  sil_out)
message("✔ Activators written: ", act_out)
message("✔ Silencers  written: ", sil_out)
message("\n✅ Peak-to-gene pipeline finished (spline correlations + raw empirical null).")
