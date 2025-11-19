# ==============================================================================
# LiftOver + Conservation scoring pipeline (rn7 -> hg38)
# Consumes the TAD-filtered LinkPeaks CSV produced earlier
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(IRanges)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(GenomeInfoDb)
  library(regioneR)
})


# ---------------------
# CONFIG — edit these
# ---------------------
CONFIG <- list(
  # Input: TAD-filtered links from the previous step
  csv_file       = "filtered_linkpeaks_all_with_TSS_noTAD.csv",
  
  # Intermediate/outputs (relative to getwd())
  prepped_csv    = "regions_for_liftover.csv",
  bed_rat_out    = "regions_rat.bed",
  
  # Chain file (rn7 -> hg38). Use the correct chain for your rn7 build.
  # Example: "C:/.../chains/rn7ToHg38.over.chain"
  chain_file     = "C:/Users/Charles/Desktop/Charles/SnMultiome/conservations/rn7ToHg38.over.chain",
  
  # LiftOver outputs
  bed_hg_out     = "regions_human.bed",
  bed_hg_merged  = "merged_regions_human.bed",
  
  # BigWig PhastCons tracks (set NA to skip)
  bw_100way      = "C:/Users/Charles/Desktop/Charles/SnMultiome/conservations/hg38.phastCons100way.bw",
  bw_7way        = "C:/Users/Charles/Desktop/Charles/SnMultiome/conservations/hg38.phastCons7way.bw",
  out_scores100  = "conservation_scores_100way.csv",
  out_scores7    = "conservation_scores_7way.csv",
  
  # Shuffled controls: provide a BED or let the script auto-generate with regioneR
  shuffled_bed   = NA,
  out_shuf100    = "shuffled_conservation_scores_100way.csv",
  out_shuf7      = "shuffled_conservation_scores_7way.csv",
  
  # Region merging (post-LiftOver, same chr; merge when gap <= merge_max_gap)
  merge_regions  = TRUE,
  merge_max_gap  = 50L,
  
  # BigWig import batching
  batch_size     = 200L
)

# ---------------------
# Helpers
# ---------------------
stop_if_missing <- function(x, msg) {
  if (is.na(x) || is.null(x) || !file.exists(x)) stop(msg, call. = FALSE)
}

ensure_ucsc_prefix <- function(chr) {
  ifelse(grepl("^chr", chr), chr, paste0("chr", chr))
}

# (0) Prepare regions from the TAD-filtered CSV
# - Ensures columns: ID, Chrx, start, end
# - Builds ID from 'peak' if present; else from coordinates
prepare_regions_from_tad_links <- function(in_csv, out_csv) {
  df <- readr::read_csv(in_csv, show_col_types = FALSE)
  
  req <- c("Chrx", "start", "end")
  if (!all(req %in% names(df))) {
    stop(sprintf("Input CSV must contain columns: %s", paste(req, collapse = ", ")), call. = FALSE)
  }
  
  df <- df |>
    mutate(
      Chrx  = ensure_ucsc_prefix(as.character(Chrx)),
      start = as.integer(start),
      end   = as.integer(end)
    )
  
  if ("peak" %in% names(df)) {
    df <- df |> mutate(ID = as.character(peak))
  } else if ("ID" %in% names(df)) {
    df <- df |> mutate(ID = as.character(ID))
  } else {
    df <- df |> mutate(ID = sprintf("pk_%s_%s_%s", Chrx, start, end))
  }
  
  bad <- which(is.finite(df$start) & is.finite(df$end) & df$end < df$start)
  if (length(bad)) {
    tmp <- df$start[bad]; df$start[bad] <- df$end[bad]; df$end[bad] <- tmp
    warning(length(bad), " rows had end<start; swapped.")
  }
  
  out <- df[, c("ID", "Chrx", "start", "end")]
  readr::write_csv(out, out_csv)
  message("[prep] Wrote: ", normalizePath(out_csv))
  out_csv
}

# Convert CSV (1-based inclusive) → GRanges (0-based half-open for BED)
as_granges_from_csv <- function(csv_file) {
  df <- readr::read_csv(csv_file, show_col_types = FALSE)
  req <- c("ID", "Chrx", "start", "end")
  if (!all(req %in% names(df))) {
    stop(sprintf("CSV must contain: %s", paste(req, collapse = ", ")), call. = FALSE)
  }
  df <- df |>
    mutate(
      Chrx  = ensure_ucsc_prefix(as.character(Chrx)),
      start = pmax(0L, as.integer(start) - 1L),  # BED 0-based start
      end   = as.integer(end),
      ID    = as.character(ID)
    )
  
  makeGRangesFromDataFrame(
    df,
    seqnames.field = "Chrx",
    start.field    = "start",
    end.field      = "end",
    keep.extra.columns = TRUE,
    starts.in.df.are.0based = TRUE
  )
}

write_bed <- function(gr, path) {
  try({ GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC" }, silent = TRUE)
  rtracklayer::export(gr, con = path, format = "BED")
}

run_liftover <- function(gr, chain_file) {
  stop_if_missing(chain_file, sprintf("Chain file not found: %s", chain_file))
  ch <- rtracklayer::import.chain(chain_file)
  
  lo <- liftOver(gr, ch)             # GRangesList
  mapped <- unlist(lo)               # All 1-to-many mapped segments
  if (length(mapped) == 0) warning("No regions mapped by LiftOver")
  
  mapped$source_index <- rep(seq_along(lo), elementNROWS(lo))
  mapped$source_ID    <- gr$ID[mapped$source_index]
  
  try({ GenomeInfoDb::seqlevelsStyle(mapped) <- "UCSC" }, silent = TRUE)
  
  cat(sprintf(
    "LiftOver: %d input → %d mapped intervals (%.1f%% with ≥1 hit)\n",
    length(gr), length(mapped), 100 * mean(elementNROWS(lo) > 0)
  ))
  mapped
}

merge_close <- function(gr, max_gap = 50L) {
  if (!length(gr)) return(gr)
  gr <- sort(gr)
  merged <- GenomicRanges::reduce(gr, min.gapwidth = max_gap + 1L)
  
  ov <- findOverlaps(merged, gr, ignore.strand = TRUE)
  merged$source_IDs <- vapply(
    split(gr$source_ID[subjectHits(ov)], queryHits(ov)),
    function(v) paste0(unique(v), collapse = ","),
    character(1)
  )
  merged$source_n <- vapply(
    split(gr$source_ID[subjectHits(ov)], queryHits(ov)),
    length,
    integer(1)
  )
  merged
}

# Weighted mean PhastCons score inside regions
weighted_means_for_batch <- function(bw_gr, regions) {
  if (!length(regions)) return(numeric())
  hits <- findOverlaps(bw_gr, regions)
  if (!length(hits)) return(rep(NA_real_, length(regions)))
  
  inter <- pintersect(bw_gr[queryHits(hits)], regions[subjectHits(hits)])
  w <- as.numeric(width(inter))
  v <- as.numeric(mcols(bw_gr)$score[queryHits(hits)])
  
  df <- data.frame(
    reg = subjectHits(hits),
    w   = w,
    v   = v
  )
  
  out <- rep(NA_real_, length(regions))
  byreg <- split(df, df$reg)
  for (k in names(byreg)) {
    d <- byreg[[k]]                # <- fix ici (]]), pas de parenthèse en trop
    out[as.integer(k)] <- sum(d$w * d$v, na.rm = TRUE) / sum(d$w, na.rm = TRUE)
  }
  out
}


extract_conservation_scores <- function(bed_path, bigwig_path, out_csv, batch_size = 200L) {
  stop_if_missing(bed_path,  sprintf("BED not found: %s", bed_path))
  stop_if_missing(bigwig_path, sprintf("BigWig not found: %s", bigwig_path))
  
  regs <- import(bed_path, format = "BED")
  try({ GenomeInfoDb::seqlevelsStyle(regs) <- "UCSC" }, silent = TRUE)
  
  n <- length(regs)
  if (n == 0) stop("No regions in BED.", call. = FALSE)
  cat(sprintf("Extracting PhastCons for %d regions from %s\n", n, basename(bigwig_path)))
  
  results <- vector("list", length = ceiling(n / batch_size))
  idx <- 1L
  for (i in seq(1L, n, by = batch_size)) {
    j <- min(i + batch_size - 1L, n)
    chunk <- regs[i:j]
    bw_gr <- import(BigWigFile(bigwig_path), which = chunk)
    means <- weighted_means_for_batch(bw_gr, chunk)
    results[[idx]] <- data.frame(
      region    = paste0(as.character(seqnames(chunk)), ":", start(chunk), "-", end(chunk)),
      mean_score= means,
      stringsAsFactors = FALSE
    )
    cat(sprintf("  batch %d: %d regions\n", idx, length(chunk)))
    idx <- idx + 1L
  }
  out <- dplyr::bind_rows(results)
  readr::write_csv(out, out_csv)
  cat("Scores written -> ", out_csv, "\n", sep = "")
  invisible(out)
}

make_shuffled <- function(gr, n = length(gr)) {
  if (!requireNamespace("regioneR", quietly = TRUE)) {
    warning("Package 'regioneR' not found; provide CONFIG$shuffled_bed to enable control.")
    return(NULL)
  }
  genome <- regioneR::getGenome("hg38")
  regioneR::randomizeRegions(gr, genome = genome, allow.overlaps = TRUE, per.chromosome = TRUE)
}

compare_distributions <- function(real_csv, shuf_csv, out_prefix) {
  real <- readr::read_csv(real_csv, show_col_types = FALSE)$mean_score
  shuf <- readr::read_csv(shuf_csv,  show_col_types = FALSE)$mean_score
  real <- real[is.finite(real)]; shuf <- shuf[is.finite(shuf)]
  
  ks <- suppressWarnings(ks.test(real, shuf))
  wt <- suppressWarnings(wilcox.test(real, shuf))
  
  df <- rbind(
    data.frame(score = real, type = "Identified"),
    data.frame(score = shuf, type = "Random")
  )
  
  # Percent histogram + dashed medians
  p <- ggplot(df, aes(score, fill = type)) +
    geom_histogram(aes(y = after_stat(count / sum(count) * 100)),
                   position = "identity", alpha = 0.5, bins = 60) +
    geom_vline(xintercept = median(real),  linetype = 2, size = 0.6) +
    geom_vline(xintercept = median(shuf),  linetype = 2, size = 0.6, color = "red") +
    scale_fill_manual(values = c(Identified = "#7aa6ff", Random = "#f08080")) +
    labs(title = "Conservation score of potential cis-regulatory regions",
         x = "Mean conservation score (%)",
         y = "Percentage of tested regions") +
    theme_bw(base_size = 12) +
    theme(legend.title = element_blank())
  
  ggsave(paste0(out_prefix, "_overlay_hist.png"), p, width = 4.5, height = 4, dpi = 300)
  
  sink(paste0(out_prefix, "_stats.txt"))
  cat("KS p=", ks$p.value, " | W p=", wt$p.value, "\n", sep = "")
  cat("median(real)=", median(real), " | median(shuf)=", median(shuf), "\n", sep = "")
  sink()
  invisible(list(ks = ks, wilcox = wt))
}

# ---------------------
# Main
# ---------------------
main <- function(cfg = CONFIG) {
  message("[0] Prepare regions from TAD-filtered LinkPeaks CSV")
  prepped <- prepare_regions_from_tad_links(cfg$csv_file, cfg$prepped_csv)
  
  message("[1] CSV -> GRanges -> BED (rat, rn7)")
  gr_rat <- as_granges_from_csv(prepped)
  write_bed(gr_rat, cfg$bed_rat_out)
  cat("BED (rn7) -> ", cfg$bed_rat_out, " (", length(gr_rat), " regions)\n", sep = "")
  
  message("[2] LiftOver rn7 -> hg38")
  gr_hg <- run_liftover(gr_rat, cfg$chain_file)
  write_bed(gr_hg, cfg$bed_hg_out)
  cat("BED (hg38) -> ", cfg$bed_hg_out, " (", length(gr_hg), " intervals)\n", sep = "")
  
  if (isTRUE(cfg$merge_regions)) {
    message("[3] Merge close hg38 regions")
    gr_merge <- merge_close(gr_hg, max_gap = cfg$merge_max_gap)
    write_bed(gr_merge, cfg$bed_hg_merged)
    cat("BED (merged) -> ", cfg$bed_hg_merged, " (", length(gr_merge), ")\n", sep = "")
  } else {
    gr_merge <- gr_hg
  }
  
  message("[4] Extract weighted PhastCons means")
  if (!is.na(cfg$bw_100way)) {
    extract_conservation_scores(cfg$bed_hg_merged, cfg$bw_100way, cfg$out_scores100, batch_size = cfg$batch_size)
  }
  if (!is.na(cfg$bw_7way)) {
    extract_conservation_scores(cfg$bed_hg_merged, cfg$bw_7way, cfg$out_scores7, batch_size = cfg$batch_size)
  }
  
  message("[5] Shuffled control (length/chr-matched)")
  if (is.na(cfg$shuffled_bed)) {
    message("  No shuffled BED provided → try regioneR to auto-generate")
    shuf <- make_shuffled(gr_merge)
    if (!is.null(shuf)) {
      cfg$shuffled_bed <- sub("\\.bed$", "_auto.bed", cfg$bed_hg_merged)
      write_bed(shuf, cfg$shuffled_bed)
      message("  Shuffled BED written: ", cfg$shuffled_bed)
    }
  }
  
  if (!is.na(cfg$shuffled_bed) && file.exists(cfg$shuffled_bed)) {
    if (!is.na(cfg$bw_100way)) {
      extract_conservation_scores(cfg$shuffled_bed, cfg$bw_100way, cfg$out_shuf100, batch_size = cfg$batch_size)
    }
    if (!is.na(cfg$bw_7way)) {
      extract_conservation_scores(cfg$shuffled_bed, cfg$bw_7way, cfg$out_shuf7, batch_size = cfg$batch_size)
    }
  }
  
  message("[6] Distribution comparisons & plots")
  if (!is.na(cfg$bw_100way) && file.exists(cfg$out_scores100) && !is.na(cfg$shuffled_bed) && file.exists(cfg$out_shuf100)) {
    compare_distributions(cfg$out_scores100, cfg$out_shuf100, out_prefix = sub("\\.csv$", "", cfg$out_scores100))
  }
  if (!is.na(cfg$bw_7way) && file.exists(cfg$out_scores7) && !is.na(cfg$shuffled_bed) && file.exists(cfg$out_shuf7)) {
    compare_distributions(cfg$out_scores7, cfg$out_shuf7, out_prefix = sub("\\.csv$", "", cfg$out_scores7))
  }
  
  message("Done.")
}

# Run when sourced via Rscript
if (sys.nframe() == 0) main()
