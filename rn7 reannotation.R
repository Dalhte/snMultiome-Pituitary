###################################################################################################################################################
#                                                                                                                                                 #
#                                                     RN7 genome reannotation                                                                     #
#                                                                                                                                                 #
###################################################################################################################################################

# --- Libraries ----------------------------------------------------------------
# Load
suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(GenomicRanges)
  library(GenomicAlignments)
  library(GenomeInfoDb)
  library(AnnotationDbi)
  library(AnnotationHub)
  library(EnsDb.Rnorvegicus.v110)
  library(BSgenome.Rnorvegicus.NCBI.rn7.2)
  library(org.Rn.eg.db)
  library(TxDb.Rnorvegicus.UCSC.rn7.refGene)
  library(rtracklayer)
  library(derfinder)
  library(GO.db)
  library(dplyr)
  library(conflicted)
  library(RMariaDB)
  library(Rsamtools)
})

# --- Paths and options --------------------------------------------------------
# Use variables for file locations so the script is easier to adapt
base_dir <- 'your/path/to/'
project_dir <- 'your/path/to/SnMultiome_Rat'
setwd(project_dir)

gtf_file <- file.path(base_dir, 'GTF', 'Rattus_norvegicus.mRatBN7.2.110.gtf')
rmsk_file <- file.path(base_dir, 'RMSK', 'rmsk.txt')
sample_bam_dir <- file.path(base_dir, 'BamSRA')


options(species = 'Rattus norvegicus')

# Source any custom functions (if available)
if (file.exists(source_function_file)) {
  source(source_function_file)
}

# --- Build TxDb from GTF -----------------------------------------------------
# Create TxDb object from GTF (Ensembl v110)
if (!file.exists(gtf_file)) stop('GTF file not found: ', gtf_file)

txdb <- makeTxDbFromGFF(gtf_file,
                        format = 'gtf',
                        dataSource = 'ensembl_V110',
                        organism = 'Rattus norvegicus')

# Obtain transcript objects and map transcript -> gene
all_tx <- transcripts(txdb, use.names = TRUE)

# Retrieve mapping from TxDb
tx_select <- AnnotationDbi::select(txdb,
                                   keys = all_tx$tx_name,
                                   keytype = 'TXNAME',
                                   columns = c('TXNAME', 'GENEID', 'TXTYPE'))

# Create named vectors for mappings
tx_to_gene <- tx_select$GENEID
names(tx_to_gene) <- tx_select$TXNAME

tx_type_map <- tx_select$TXTYPE
names(tx_type_map) <- tx_select$TXNAME

# Attach gene ids to transcript metadata where available
mcols(all_tx)$gene_id <- tx_to_gene[mcols(all_tx)$tx_name]

# Transcript ends (1-bp ranges at transcript end)
tx_ends <- resize(all_tx, width = 1, fix = 'end')

# --- Extract basic genomic features -----------------------------------------
# 3' UTRs, 5' UTRs, exons and introns by transcript (unlisted for convenience)
three_utrs <- unlist(threeUTRsByTranscript(txdb, use.names = TRUE), use.names = TRUE)
mcols(three_utrs)$gene_id <- tx_to_gene[names(three_utrs)]

five_utrs <- unlist(fiveUTRsByTranscript(txdb, use.names = TRUE), use.names = TRUE)
mcols(five_utrs)$gene_id <- tx_to_gene[names(five_utrs)]

exons <- unlist(exonsBy(txdb, by = 'tx', use.names = TRUE), use.names = TRUE)
mcols(exons)$gene_id <- tx_to_gene[names(exons)]

introns <- unlist(intronsByTranscript(txdb, use.names = TRUE), use.names = TRUE)
mcols(introns)$gene_id <- tx_to_gene[names(introns)]

# Export BED files for features (optional; kept for reproducibility)
rtracklayer::export.bed(three_utrs,  './threeUTRs.bed')
rtracklayer::export.bed(five_utrs,   './fiveUTRs.bed')
rtracklayer::export.bed(exons,       './exons.bed')
rtracklayer::export.bed(introns,     './introns.bed')

# --- Coverage and region detection with derfinder ----------------------------
# Identify BAM files for plus and minus strands using a filename pattern
get_bam_files <- function(dir, pattern) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  return(sort(files))
}

files_plus <- get_bam_files(sample_bam_dir, 'MergedRN7.1fwd.bam')
files_minus <- get_bam_files(sample_bam_dir, 'MergedRN7.1rev.bam')

if (length(files_plus) == 0 && length(files_minus) == 0) {
  warning('No BAM files found with the given patterns in: ', sample_bam_dir)
}

# Chromosomes of interest (NCBI-style names without 'chr' prefix)
chrs <- c(as.character(1:20), 'X', 'Y')

# Compute full coverage objects (cutoff chosen as 25 in original script)
full_cov_plus <- fullCoverage(files = files_plus, chrs = chrs, verbose = FALSE,
                              txdb = txdb, cutoff = 25, chrsStyle = 'NCBI')

full_cov_minus <- fullCoverage(files = files_minus, chrs = chrs, verbose = FALSE,
                               txdb = txdb, cutoff = 25, chrsStyle = 'NCBI')

# Total mapped reads per sample (vector)
total_mapped_plus <- sapply(files_plus, getTotalMapped, chrs = chrs)
total_mapped_minus <- sapply(files_minus, getTotalMapped, chrs = chrs)

# Build region matrices. Keep original parameters (targetSize and maxRegionGap)
region_mat_plus <- regionMatrix(full_cov_plus, cutoff = 25, L = 1,
                                maxRegionGap = 50L, totalMapped = total_mapped_plus,
                                targetSize = 31000000, txdb = txdb, chrsStyle = 'NCBI')

region_mat_minus <- regionMatrix(full_cov_minus, cutoff = 25, L = 1,
                                 maxRegionGap = 50L, totalMapped = total_mapped_minus,
                                 targetSize = 31000000, txdb = txdb, chrsStyle = 'NCBI')

# Extract GRanges lists for detected regions and set strand
region_list_plus <- lapply(region_mat_plus, '[[', 1)
region_pos <- unlist(GRangesList(region_list_plus), use.names = TRUE)
strand(region_pos) <- '+'

region_list_minus <- lapply(region_mat_minus, '[[', 1)
region_neg <- unlist(GRangesList(region_list_minus), use.names = TRUE)
strand(region_neg) <- '-'

# Combine and slightly extend detected regions upstream (50 bp)
regions_all_original <- c(region_pos, region_neg)
regions_extended <- grange_extend(regions_all_original, upstream = 50, downstream = 0)

# Merge regions that overlap or are adjacent
regions_merged <- reduce(regions_extended)
names(regions_merged) <- as.character(seq_along(regions_merged))
regions_merged$region_name <- names(regions_merged)

# Export merged regions
rtracklayer::export.bed(regions_merged, './regions_all.bed')

# --- Annotate regions against gene features ---------------------------------
# Build mapping of regions to feature lists (3' UTRs, exons, 5' UTRs, introns)
feature_list <- list(threeUTRs = three_utrs,
                     exons     = exons,
                     fiveUTRs  = five_utrs,
                     introns   = introns)

map_to_features <- lapply(feature_list, function(feat) unique(from(findOverlaps(regions_merged, feat))))
map_to_features$intergenic <- setdiff(seq_along(regions_merged), unique(unlist(map_to_features)))

# Quick summary table: binary matrix of annotations per region
regions_annotation_matrix <- as.data.frame(do.call(cbind, lapply(map_to_features, function(idx) ifelse(seq_along(regions_merged) %in% idx, 1, 0))))

# --- RepeatMasker overlap ----------------------------------------------------
# Read RepeatMasker text table and convert to GRanges
if (file.exists(rmsk_file)) {
  rmsk_txt <- read.csv(rmsk_file, sep = '\t', header = FALSE, stringsAsFactors = FALSE)
  rmsk_df <- data.frame(
    seqnames = gsub('.*chr', '', rmsk_txt$V6),
    start    = rmsk_txt$V7,
    end      = rmsk_txt$V8,
    strand   = rmsk_txt$V10,
    gene_name = rmsk_txt$V11,
    stringsAsFactors = FALSE
  )
  rmsk_gr <- makeGRangesFromDataFrame(rmsk_df, keep.extra.columns = TRUE)
  rtracklayer::export.bed(rmsk_gr, './rmskRN7.bed')
  
  map_to_rmsk <- unique(from(findOverlaps(regions_merged, rmsk_gr, select = 'all')))
  # Consider regions within 20 bp downstream of repeat ends as associated
  rmsk_ends <- resize(rmsk_gr, width = 1, fix = 'end')
  nearest_to_rmsk <- distanceToNearest(regions_merged, rmsk_ends, select = 'arbitrary')
  close_to_rmsk <- from(subset(nearest_to_rmsk, mcols(nearest_to_rmsk)$distance < 20))
  map_to_rmsk <- unique(c(map_to_rmsk, close_to_rmsk))
} else {
  warning('RepeatMasker file not found at: ', rmsk_file)
  map_to_rmsk <- integer(0)
}

# --- Map regions to transcripts and annotate intronic/intergenic signals ----
# Map regions labelled as exonic to transcripts and extract gene ids
exonic_idx <- map_to_features$exons
if (length(exonic_idx) > 0) {
  hits_exons <- findOverlaps(regions_merged[exonic_idx], all_tx)
  anno_exons <- reformat_anno(hits_exons, all_tx, type = 'tx', regions_merged[exonic_idx], distance = NA)
  exons_genes <- unique(as.character(anno_exons$gene_id))
} else {
  anno_exons <- NULL
  exons_genes <- character(0)
}

# Intronic regions (excluding those overlapping 3'UTR, exons, 5'UTR or repeats)
intronic_idx <- setdiff(map_to_features$introns, c(map_to_features$threeUTRs, map_to_features$exons, map_to_features$fiveUTRs, map_to_rmsk))
if (length(intronic_idx) > 0) {
  intronic_regions <- regions_merged[intronic_idx]
  hits_intronic <- findOverlaps(intronic_regions, all_tx, select = 'all')
  dist_to_exon <- distanceToNearest(intronic_regions[from(hits_intronic)], exons, select = 'arbitrary')
  anno_intronic <- reformat_anno(hits_intronic, all_tx, type = 'tx', intronic_regions, distance = mcols(dist_to_exon)$distance)
  intronic_genes <- unique(anno_intronic$gene_id)
} else {
  anno_intronic <- NULL
  intronic_genes <- character(0)
}

# Intergenic regions (excluding repeats)
intergenic_idx <- setdiff(map_to_features$intergenic, map_to_rmsk)
intergenic_regions <- regions_merged[intergenic_idx]

# For intergenic regions, find nearest transcript ends (respecting strand)
if (length(intergenic_regions) > 0) {
  intergenic_to_ends <- follow(intergenic_regions, tx_ends, select = 'all', ignore.strand = FALSE)
  dist_to_ends <- distance(intergenic_regions[queryHits(intergenic_to_ends)], tx_ends[subjectHits(intergenic_to_ends)])
  anno_ends <- reformat_anno(intergenic_to_ends, tx_ends, type = 'ends', intergenic_regions, distance = dist_to_ends)
  # Keep plausible downstream annotations: transcript type == 'transcript' and distance < 1000bp
  mcols(anno_ends)$geneLength <- width(genes(txdb)) [mcols(anno_ends)$gene_id]
  mcols(anno_ends)$genetype <- tx_type_map[mcols(anno_ends)$tx_name]
  anno_ends$perc <- with(mcols(anno_ends), distance / geneLength)
  anno_ends_valid <- subset(anno_ends, genetype == 'transcript' & distance < 1000)
  ends_genes <- unique(anno_ends_valid$gene_id)
} else {
  anno_ends_valid <- NULL
  ends_genes <- character(0)
}

# --- Collect novel candidate regions and export --------------------------------
# Build GRanges for novel intronic, exonic and downstream (ends) signals
# Prepare intronic novel regions
if (!is.null(anno_intronic) && nrow(anno_intronic) > 0) {
  intronic_df <- anno_intronic[, c('region_name', 'tx_name', 'gene_id', 'type')]
  intronic_gr <- makeGRangesFromDataFrame(intronic_df, keep.extra.columns = TRUE)
  intronic_novel <- unique(grange_extend(intronic_gr, -45, 5)[, c('tx_name', 'gene_id')])
  mcols(intronic_novel)$gene_name <- tx_to_gene[intronic_novel$tx_name]
} else {
  intronic_novel <- GRanges()
}

# Prepare downstream novel ends
if (!is.null(anno_ends_valid) && nrow(anno_ends_valid) > 0) {
  ends_gr <- makeGRangesFromDataFrame(anno_ends_valid, keep.extra.columns = TRUE)
  ends_resized <- resize(ends_gr, width = 1, fix = 'end')
  upstream_txend_hits <- follow(ends_resized, tx_ends, select = 'last', ignore.strand = FALSE)
  ends_df <- cbind(as.data.frame(ends_resized),
                   txend = start(tx_ends[upstream_txend_hits, drop = TRUE]))
  # Create a combined GRanges spanning txend and detected position
  combined_df <- data.frame(seqnames = ends_df$seqnames,
                            start    = pmin(as.integer(ends_df$start), as.integer(ends_df$txend)),
                            end      = pmax(as.integer(ends_df$end),   as.integer(ends_df$txend)),
                            strand   = ends_df$strand,
                            gene_id  = ends_df$gene_id,
                            tx_name  = ends_df$tx_name,
                            gene_name = ends_df$gene_name,
                            stringsAsFactors = FALSE)
  ends_novel <- makeGRangesFromDataFrame(combined_df, keep.extra.columns = TRUE)
} else {
  ends_novel <- GRanges()
}

# Prepare exonic novel regions
if (!is.null(anno_exons) && nrow(anno_exons) > 0) {
  exon_df <- as.data.frame(anno_exons)
  exon_df_keep <- data.frame(seqnames = exon_df$seqnames,
                             start = exon_df$start,
                             end   = exon_df$end,
                             strand = exon_df$strand,
                             gene_id = exon_df$gene_id,
                             tx_name = exon_df$tx_name,
                             type = exon_df$type,
                             stringsAsFactors = FALSE)
  exon_novel <- makeGRangesFromDataFrame(exon_df_keep, keep.extra.columns = TRUE)
} else {
  exon_novel <- GRanges()
}

# Tag types and combine novel UTR candidates
if (length(intronic_novel) > 0) mcols(intronic_novel)$type <- 'protein_coding'
if (length(ends_novel) > 0) mcols(ends_novel)$type <- 'protein_coding'
if (length(exon_novel) > 0) mcols(exon_novel)$type <- 'protein_coding'

novel_utrs_all <- c(exon_novel, intronic_novel, ends_novel)
if (length(novel_utrs_all) > 0) mcols(novel_utrs_all)$transcript_id <- novel_utrs_all$tx_name

# Add gene names from original GTF where possible
if (exists('essai2')) rm(essai2) # avoid accidental reuse
essai2 <- rtracklayer::import(gtf_file)
name_map <- setNames(essai2$gene_name, essai2$gene_id)
if (length(novel_utrs_all) > 0) {
  mcols(novel_utrs_all)$gene_name <- name_map[novel_utrs_all$gene_id]
}

# Export novel candidates as BED and create a GTF/GFF-like output table
if (length(novel_utrs_all) > 0) {
  rtracklayer::export.bed(novel_utrs_all, './novel_utrs_all.bed')
  
  df_novel <- as.data.frame(novel_utrs_all, row.names = NULL)
  df_novel$uniq_id <- seq_len(nrow(df_novel))
  df_novel$transcript_name <- paste(df_novel$gene_name, df_novel$uniq_id, sep = '.')
  
  gtf_out <- data.frame(
    chr = as.character(seqnames(novel_utrs_all)),
    source = 'DL',
    feature = 'exon',
    start = start(novel_utrs_all),
    end = end(novel_utrs_all),
    score = '.',
    strand = as.character(strand(novel_utrs_all)),
    frame = '.',
    attribute = paste0('gene_id "', df_novel$gene_id, '"; ',
                       'gene_name "', df_novel$gene_name, '"; ',
                       'transcript_id "', df_novel$transcript_id, '"; ',
                       'transcript_name "', df_novel$transcript_name, '"; ',
                       'gene_biotype "protein_coding"; ',
                       'transcript_biotype "protein_coding";')
  )
  
  write.table(unique(gtf_out), './NovelPitAnnotLevelRN7.gff3', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
}

# --- Summary lists -----------------------------------------------------------
exonic_genes_list   <- if (exists('exons_genes')) exons_genes else character(0)
intronic_genes_list <- if (exists('intronic_genes')) intronic_genes else character(0)
downstream_genes    <- if (exists('ends_genes')) ends_genes else character(0)

datalist <- list(
  'Exonic signal'      = exonic_genes_list,
  'Intronic signal'    = intronic_genes_list,
  'Downstream signal'  = downstream_genes
)

