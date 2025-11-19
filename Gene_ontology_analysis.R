# =============================================================================
# Script 11 — Pathway enrichment (GO-ALL / KEGG / Hallmark / Reactome)
#              and ontology dotplot per expression cluster
# =============================================================================
# Requirements:
#   - RDS file from Script 10: segmentation_results.rds (contains res_expr)
#   - Expression matrix (genes × pseudotime bins)
# Output:
#   - results/pseudotime_enrichment/enrichments_*.csv
#   - results/pseudotime_enrichment/dotplot_by_cluster.svg
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
  library(msigdbr)
  library(AnnotationDbi)
  library(org.Rn.eg.db)
  library(GOSemSim)
  library(tidytext)
  library(text2vec)
  library(purrr)
})

# -------------------------- User parameters ----------------------------------
expr_matrix_path <- "outputs/regtime/data_regtime.csv"          # gene × bin
res_rds_path     <- "outputs/heatmaps/segmentation_results.rds" # from Script 10
out_dir          <- "results/pseudotime_enrichment"
gene_pattern_exclude <- c("^LOC","^NEWGENE","^RGD")
min_gs_size <- 5
q_cut       <- 0.10
set.seed(42)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------- Load res_expr ------------------------------------
if (!exists("res_expr")) {
  if (file.exists(res_rds_path)) {
    seg_res <- readRDS(res_rds_path)
    if (is.list(seg_res) && "res_expr" %in% names(seg_res)) {
      res_expr <- seg_res$res_expr
    } else if (is.list(seg_res) && !is.null(seg_res$mat_ord)) {
      res_expr <- seg_res
    } else {
      stop("res_expr object not found inside ", res_rds_path)
    }
  } else stop("res_expr not loaded and no RDS found at ", res_rds_path)
}
stopifnot(is.list(res_expr), !is.null(res_expr$mat_ord), !is.null(res_expr$seg_id))

# -------------------------- Universe (expression matrix) ---------------------
expr_mat <- readr::read_csv(expr_matrix_path, col_types = readr::cols()) %>%
  dplyr::rename(Gene = 1) %>%
  tibble::column_to_rownames("Gene") %>%
  as.matrix()
expr_mat <- expr_mat[!grepl(paste(gene_pattern_exclude, collapse="|"),
                            rownames(expr_mat)), , drop = FALSE]
expr_mat[is.na(expr_mat)] <- 0

expr_universe <- bitr(rownames(expr_mat), fromType="SYMBOL", toType="ENTREZID",
                      OrgDb = org.Rn.eg.db) %>% pull(ENTREZID) %>% unique()

# -------------------------- Expression clusters ------------------------------
rownames_ok <- intersect(rownames(res_expr$mat_ord), rownames(expr_mat))
stopifnot(length(rownames_ok) > 0)

expr_gene2seg <- tibble(
  SYMBOL  = rownames(res_expr$mat_ord),
  Cluster = paste0("E", res_expr$seg_id)
) %>% filter(SYMBOL %in% rownames_ok)

genes_by_cluster <- split(expr_gene2seg$SYMBOL, expr_gene2seg$Cluster)

# Convert SYMBOL → ENTREZ
genes_list_rat <- lapply(genes_by_cluster, function(symb) {
  out <- bitr(symb, fromType="SYMBOL", toType="ENTREZID", OrgDb = org.Rn.eg.db)
  unique(out$ENTREZID)
})

# -------------------------- Over-representation analysis ---------------------
run_cc <- function(genes, fun, label, ...) {
  compareCluster(
    geneCluster   = genes,
    fun           = fun,
    minGSSize     = min_gs_size,
    maxGSSize     = 500,
    pAdjustMethod = "BH",
    qvalueCutoff  = q_cut,
    ...
  )
}

cmp_GO_ALL <- run_cc(
  genes_list_rat, "enrichGO", "GO (rat)",
  OrgDb = org.Rn.eg.db, keyType = "ENTREZID", ont = "ALL",
  universe = expr_universe
)
cmp_KEGG <- run_cc(
  genes_list_rat, "enrichKEGG", "KEGG (rat)",
  organism = "rno", keyType = "ncbi-geneid",
  universe = expr_universe
)
hall_tbl <- msigdbr(species = "Rattus norvegicus", category = "H")
hall_T2G <- hall_tbl %>% dplyr::select(gs_name, ncbi_gene) %>% distinct()
cmp_HALLMARK <- run_cc(
  genes_list_rat, "enricher", "Hallmark (rat)",
  TERM2GENE = hall_T2G, universe = expr_universe
)
react_tbl <- msigdbr(species = "Rattus norvegicus", category = "C2", subcategory = "CP:REACTOME")
react_T2G <- react_tbl %>% dplyr::select(gs_name, ncbi_gene) %>% distinct()
cmp_REACT <- run_cc(
  genes_list_rat, "enricher", "Reactome (rat)",
  TERM2GENE = react_T2G, universe = expr_universe
)

# -------------------------- Aggregate all enrichment results -----------------
bind_safe <- function(cc, label) {
  if (is.null(cc)) return(NULL)
  if (nrow(cc@compareClusterResult) == 0) return(NULL)
  cc@compareClusterResult %>% mutate(Source = label)
}
enrich_tbl <- list(
  bind_safe(cmp_GO_ALL,   "GO (rat)"),
  bind_safe(cmp_KEGG,     "KEGG (rat)"),
  bind_safe(cmp_HALLMARK, "Hallmark (rat)"),
  bind_safe(cmp_REACT,    "Reactome (rat)")
) %>% bind_rows()
stopifnot(nrow(enrich_tbl) > 0)

# -------------------------- Compact redundancy removal -----------------------
options(Matrix.warnDeprecatedCoerce = FALSE)

make_er <- function(df, ontology = "BP") {
  new("enrichResult",
      result       = distinct(df, ID, .keep_all = TRUE),
      ontology     = ontology,
      pvalueCutoff = 0.05, qvalueCutoff = 0.20,
      organism     = "rat", keytype = "ENTREZID")
}

tmp_split <- enrich_tbl %>%
  mutate(Cluster = as.character(Cluster)) %>%
  group_by(Cluster, Source) %>%
  group_split()
nm <- paste0(
  map_chr(tmp_split, ~ unique(.x$Cluster)[1]), ".",
  map_chr(tmp_split, ~ unique(.x$Source)[1])
)
enrich_ls <- map(tmp_split, ~ make_er(.x, ontology = unique(.x$ONTOLOGY %||% "BP")))
names(enrich_ls) <- nm

cluster2genes_tbl <- enrich_tbl %>%
  dplyr::select(Cluster, geneID) %>%
  separate_rows(geneID, sep = "/") %>%
  mutate(geneID = str_trim(geneID)) %>%
  filter(geneID != "") %>% distinct()
cluster2genes <- split(cluster2genes_tbl$geneID, cluster2genes_tbl$Cluster) %>% lapply(unique)

get_sem <- local({
  cache <- list()
  function(ont) {
    ont <- as.character(ont)[1]
    if (is.na(ont) || !ont %in% c("BP", "CC", "MF")) return(NULL)
    if (!ont %in% names(cache))
      cache[[ont]] <<- tryCatch(godata("org.Rn.eg.db", ont, computeIC = TRUE), error = function(e) NULL)
    cache[[ont]]
  }
})

safe_z <- function(k, M, n, N) {
  if (anyNA(c(k, M, n, N)) || any(c(k, M, n, N) < 0) || (M + N) <= 1) return(0)
  den <- n * M * N * (M + N - n)
  den <- den / ((M + N)^2 * (M + N - 1))
  if (!is.finite(den) || den == 0) return(0)
  z <- (k - n * M / (M + N)) / sqrt(den)
  if (!is.finite(z)) 0 else z
}

fuse_block <- function(df_block, cl_genes) {
  # select one representative term per block (break ties deterministically)
  best <- dplyr::slice_min(df_block, p.adjust, n = 1, with_ties = FALSE)
  
  # parse background and sample sizes from that single row
  k <- as.integer(best$Count[[1]])
  M <- suppressWarnings(as.integer(stringr::str_extract(best$BgRatio[[1]], "^[0-9]+")))
  U <- suppressWarnings(as.integer(stringr::str_extract(best$BgRatio[[1]], "(?<=/)[0-9]+$")))
  n <- suppressWarnings(as.integer(stringr::str_extract(best$GeneRatio[[1]], "(?<=/)[0-9]+$")))
  
  # fallbacks if parsing failed
  if (!is.finite(M) || is.na(M)) M <- 0L
  if (!is.finite(U) || is.na(U)) U <- 0L
  if (!is.finite(n) || is.na(n)) n <- length(cl_genes)
  
  # one-tailed hypergeometric and derived metrics
  p <- stats::phyper(k - 1, M, U - M, n, lower.tail = FALSE)
  
  dplyr::mutate(
    best,
    pvalue         = p,
    p.adjust       = p.adjust(p, method = "BH"),
    qvalue         = p.adjust,
    RichFactor     = ifelse(M > 0, k / M, NA_real_),
    FoldEnrichment = ifelse(M > 0 && n > 0 && U > 0, (k / n) / (M / U), NA_real_),
    zScore         = safe_z(k, M, n, max(U - M, 0L))
  )
}
collapse_one <- function(er, wang_cut = 0.6, jac_cut = 0.4, lv_dist = 6) {
  # Requires: cluster2genes (list), get_sem(), fuse_block() in scope.
  
  cl <- unique(er@result$Cluster)[1]
  cl_genes <- cluster2genes[[cl]]
  if (is.null(cl_genes)) cl_genes <- character(0)
  
  # -- (a) Semantic reduction (Wang) + simplify (GO only) --------------------
  sem <- get_sem(er@ontology)
  er1 <- if (!is.null(sem) && nrow(er@result) > 1) {
    er_tmp <- tryCatch(
      clusterProfiler::pairwise_termsim(er, method = "Wang", semData = sem),
      error = function(e) er
    )
    clusterProfiler::simplify(er_tmp, cutoff = wang_cut, by = "p.adjust", select_fun = min)
  } else {
    er
  }
  
  # Helper: sanitize a similarity matrix into a valid distance + hclust
  build_hc <- function(S) {
    S <- as.matrix(S)
    if (!is.matrix(S) || nrow(S) < 2) return(NULL)
    S[!is.finite(S)] <- 0
    # force symmetry + unit diagonal and clamp to [0,1]
    S <- (S + t(S)) / 2
    diag(S) <- 1
    S <- pmin(pmax(S, 0), 1)
    D <- 1 - S
    diag(D) <- 0
    if (any(!is.finite(D))) return(NULL)
    if (sum(D, na.rm = TRUE) <= 0) return(NULL)
    hc <- stats::hclust(stats::as.dist(D), method = "average")
    # enforce non-decreasing heights (fixes 'height not sorted' error)
    hc$height <- cummax(hc$height + .Machine$double.eps)
    hc
  }
  
  # -- (b) Jaccard over gene sets -------------------------------------------
  er2 <- er1
  if (nrow(er1@result) > 1) {
    er2 <- tryCatch(
      suppressWarnings(clusterProfiler::pairwise_termsim(er1, method = "jc")),
      error = function(e) er1
    )
    if (!is.null(er2@termsim)) {
      # normalize sparse storage
      if (inherits(er2@termsim, "dgTMatrix")) {
        er2@termsim <- methods::as(er2@termsim, "CsparseMatrix")
      }
      hc <- build_hc(er2@termsim)
      if (!is.null(hc)) {
        grp <- stats::cutree(hc, h = jac_cut)
        # map group labels back to result rows
        names(grp) <- rownames(as.matrix(er2@termsim))
        er2@result <- er2@result |>
          dplyr::mutate(grp = grp[Description]) |>
          dplyr::group_split(grp) |>
          purrr::map_dfr(fuse_block, cl_genes = cl_genes)
      }
    }
  }
  
  # -- (c) Text Jaccard on term descriptions --------------------------------
  er3 <- er2
  if (nrow(er2@result) > 1) {
    bow <- er2@result |>
      dplyr::mutate(id = dplyr::row_number()) |>
      tidytext::unnest_tokens(word, Description) |>
      dplyr::anti_join(tidytext::stop_words, by = "word") |>
      dplyr::filter(!stringr::str_detect(word, "^[0-9]+$"))
    if (nrow(bow) > 1) {
      dtm <- bow |>
        dplyr::count(id, word, name = "n") |>
        tidytext::cast_sparse(id, word, n)
      dtm <- methods::as(dtm, "CsparseMatrix")
      
      S <- suppressWarnings(text2vec::sim2(dtm, dtm, method = "jaccard", norm = "none"))
      hc <- build_hc(S)
      if (!is.null(hc)) {
        grp <- stats::cutree(hc, h = lv_dist / 10)
        er3@result <- er2@result |>
          dplyr::mutate(grp = grp) |>
          dplyr::group_split(grp) |>
          purrr::map_dfr(fuse_block, cl_genes = cl_genes)
      }
    }
  }
  
  er3
}


collapsed_ls <- imap(enrich_ls, ~ collapse_one(.x, .6, .4, 6))
collapsed_tbl <- imap_dfr(collapsed_ls, ~ .x@result %>% mutate(ClusterSource = .y)) %>%
  arrange(Cluster, p.adjust)

# -------------------------- Recompute consistent statistics ------------------
enrich_sizes <- enrich_tbl %>%
  transmute(
    Cluster, ID,
    M_orig = as.integer(str_extract(BgRatio, "^[0-9]+")),
    U_orig = as.integer(str_extract(BgRatio, "(?<=/)[0-9]+$"))
  ) %>% distinct()

collapsed_tbl_fixed <- collapsed_tbl %>%
  left_join(enrich_sizes, by = c("Cluster", "ID"), relationship = "many-to-many") %>%
  rowwise() %>%
  mutate(
    k = Count,
    n = as.integer(str_extract(GeneRatio, "(?<=/)[0-9]+$")),
    M = M_orig, U = U_orig,
    GeneRatio      = sprintf("%d/%d", k, n),
    BgRatio        = sprintf("%d/%d", M, U),
    RichFactor     = k / M,
    FoldEnrichment = (k / n) / (M / U),
    zScore         = safe_z(k, M, n, U - M),
    pvalue         = phyper(k - 1, M, U - M, n, lower.tail = FALSE)
  ) %>% ungroup() %>%
  group_by(Cluster) %>%
  mutate(p.adjust = p.adjust(pvalue, method = "BH"), qvalue = p.adjust) %>%
  ungroup() %>%
  dplyr::select(-M_orig, -U_orig, -M, -U, -k, -n)

# Save tables
write_csv(enrich_tbl,         file.path(out_dir, "enrichments_by_cluster.csv"))
write_csv(collapsed_tbl,      file.path(out_dir, "enrichments_collapsed_union.csv"))
write_csv(collapsed_tbl_fixed,file.path(out_dir, "enrichments_collapsed_union_FIXED.csv"))

# -------------------------- Dotplot generation -------------------------------
df0 <- collapsed_tbl_fixed %>%
  mutate(
    label = str_to_sentence(Description) %>% str_wrap(60),
    n   = as.numeric(str_extract(GeneRatio, "(?<=/)[0-9]+")),
    GeneRatio_num = Count / n,
    score = -log10(qvalue + 1e-8)
  )

df_uni <- df0 %>%
  group_by(Cluster, label) %>%
  slice_max(FoldEnrichment, n = 1, with_ties = FALSE) %>%
  ungroup()

top_labels <- df_uni %>%
  group_by(Cluster) %>%
  slice_max(FoldEnrichment, n = 20, with_ties = FALSE) %>%
  pull(label) %>% unique()

plot_df <- df_uni %>% filter(label %in% top_labels)

order_vec <- character(0)
for (cl in unique(plot_df$Cluster)) {
  new_lab <- plot_df %>%
    filter(Cluster == cl) %>%
    arrange(desc(FoldEnrichment)) %>%
    pull(label)
  order_vec <- c(order_vec, setdiff(new_lab, order_vec))
}
plot_df <- plot_df %>%
  mutate(label_ord = factor(label, levels = rev(order_vec)))

max_GR <- max(plot_df$GeneRatio_num, na.rm = TRUE)
size_breaks <- c(0.05, 0.10, 0.15, 0.20, 0.3, 0.4)
size_breaks <- size_breaks[size_breaks <= max_GR]
if (length(size_breaks) < 3) {
  size_breaks <- scales::pretty_breaks(n = 5)(c(0, max_GR))
  size_breaks <- size_breaks[size_breaks > 0 & size_breaks <= max_GR]
}

p_dot <- ggplot(plot_df, aes(Cluster, label_ord, size = GeneRatio_num, fill = score)) +
  geom_point(shape = 21, colour = "black", stroke = .25, alpha = .9) +
  scale_fill_gradientn(colours = colorRampPalette(c("blue", "red"))(100),
                       name = expression(-log[10](q))) +
  scale_size(range = c(2, 7), breaks = size_breaks,
             labels = scales::percent_format(accuracy = 0.1), name = "GeneRatio") +
  labs(title = "Top enriched biological processes per cluster",
       subtitle = "One term per Cluster × label (highest fold-enrichment)",
       x = NULL, y = NULL) +
  theme_bw(base_size = 10) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 7))

dot_svg <- file.path(out_dir, "dotplot_by_cluster.svg")
ggsave(dot_svg, p_dot, width = 6, height = 12, units = "in")
message("Dotplot saved: ", dot_svg)
