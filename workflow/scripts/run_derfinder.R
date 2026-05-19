# Coverage-based expressed-region caller, used as a baseline against fastder.
#
# Reads every per-sample BigWig in --bigwig-dir, normalises each sample's
# coverage to CPM (counts per million) using the same library_size formula
# fastder applies internally, averages the per-sample CPMs across all
# samples (zero-included for samples without coverage at a position), and
# emits one GTF transcript per maximal run of consecutive bases whose mean
# CPM is at or above --cutoff.
#
# The semantics match what `derfinder::findRegions` does on a single
# coverage track, but we extract IRanges from the Rle directly so the
# script is robust to derfinder API changes (older versions of findRegions
# choke on logical Rles built from coverage that contains NA bases, even
# after scrubbing).
#
# Output schema (per region): one gene + one transcript + one exon line,
# each scored by the mean CPM in that region.
#
# Parameter equivalence with fastder:
#   --cutoff       <-> fastder --min-coverage      (CPM)
#   --min-length   <-> fastder --min-length        (bp; post-filter)
#   --maxregiongap <-> fastder --position-tolerance (bp; gap-bridging slack)
#   --chromosomes  <-> fastder --chr               (used for library_size scope)
#
# The --maxregiongap mapping is a behavioural analogue rather than an
# exact equivalent: derfinder bridges below-threshold gaps inside a
# region, while fastder allows a few bp of slack at SJ-anchored
# boundaries. Both are tolerance knobs in the bp dimension; the report
# sweeps both to surface where they end up roughly comparable.
suppressPackageStartupMessages({
  library(optparse)
  library(derfinder)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(rtracklayer)
})

opt_list <- list(
  make_option("--bigwig-dir", type = "character", help = "Directory holding per-sample .bw files"),
  make_option("--out-gtf", type = "character", help = "Path of the GTF to write"),
  make_option("--cutoff", type = "double", default = 0.05,
              help = "Coverage cutoff above which a base counts as expressed (in CPM)"),
  make_option("--chromosomes", type = "character", default = NULL,
              help = "Space-separated list of chromosomes to analyse; default = all in the BigWigs"),
  make_option("--maxregiongap", type = "integer", default = 0,
              help = "Bridge below-threshold gaps of at most this many bases inside a region"),
  make_option("--min-length", type = "integer", default = 1,
              help = "Drop regions shorter than this many bases")
)
parser <- OptionParser(option_list = opt_list)
arg_strings <- commandArgs(trailingOnly = TRUE)
chrom_idx <- match("--chromosomes", arg_strings)
chroms <- character(0)
if (!is.na(chrom_idx)) {
  i <- chrom_idx + 1
  while (i <= length(arg_strings) && !startsWith(arg_strings[i], "--")) {
    chroms <- c(chroms, arg_strings[i])
    i <- i + 1
  }
  arg_strings <- arg_strings[-(chrom_idx:(i - 1))]
}
opt <- parse_args(parser, args = arg_strings)
if (length(chroms) == 0 && !is.null(opt$chromosomes)) {
  chroms <- strsplit(opt$chromosomes, "[ ,]+")[[1]]
}

bw_files <- list.files(opt$`bigwig-dir`, pattern = "\\.(all|plus|minus)\\.bw$",
                       full.names = TRUE)
if (length(bw_files) == 0) {
  stop("[run_derfinder] No .bw files in ", opt$`bigwig-dir`)
}
# A stranded sample has a .plus.bw and a .minus.bw; both strand tracks
# belong to one sample. Group them so the two strands count as one sample,
# not two, in the library-size and averaging steps below.
sample_of <- sub("\\.(all|plus|minus)\\.bw$", "", basename(bw_files))
samples <- split(unname(bw_files), sample_of)
sample_ids <- names(samples)
flat_files <- unlist(samples, use.names = FALSE)
file_sample_idx <- rep(seq_along(samples), lengths(samples))

if (length(chroms) == 0) {
  chroms <- seqlevels(rtracklayer::BigWigFile(flat_files[1]))
}

# library_size per sample, scoped to the user's chromosomes. Matches
# fastder's accumulator: total_reads = (end - start) * value, summed over
# only the chromosomes passed via --chr.
compute_library_size <- function(bw_path, target_chroms) {
  available <- seqlevels(rtracklayer::BigWigFile(bw_path))
  keep <- intersect(target_chroms, available)
  if (length(keep) == 0) return(0)
  gr <- rtracklayer::import.bw(bw_path,
                               which = GRanges(seqnames = keep,
                                               ranges = IRanges(1L, .Machine$integer.max %/% 2L)))
  if (length(gr) == 0) return(0)
  sum(as.numeric(width(gr)) * as.numeric(gr$score))
}

message("[run_derfinder] ", length(flat_files), " BigWig files grouped into ",
        length(samples), " samples; computing library sizes")
lib_sizes <- vapply(samples, function(files)
  sum(vapply(files, compute_library_size, numeric(1), target_chroms = chroms)),
  FUN.VALUE = numeric(1))
cpm_factors <- lib_sizes / 1e6
for (i in seq_along(samples)) {
  if (cpm_factors[i] <= 0) {
    message("[run_derfinder] WARN: sample ", sample_ids[i],
            " has empty library_size on chromosomes ",
            paste(chroms, collapse = ","), "; sample will be skipped.")
  }
}

regions_per_chrom <- list()
for (chrom in chroms) {
  message("[run_derfinder] processing ", chrom)
  fullCov <- tryCatch(
    fullCoverage(files = flat_files, chrs = chrom, verbose = FALSE),
    error = function(e) {
      message("[run_derfinder]   skip ", chrom, ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(fullCov) || length(fullCov[[chrom]]) == 0) next

  # Mean per-sample CPM across samples. fullCoverage returns one Rle per
  # BigWig file; a stranded sample's plus and minus tracks are summed back
  # together here so each sample counts once. Bases unobserved by all
  # BigWigs come back as NA; treat them as zero (NA in the threshold mask
  # would propagate into IRanges construction).
  per_file_rles <- as.list(fullCov[[chrom]])
  acc <- numeric(0)
  for (si in seq_along(samples)) {
    if (cpm_factors[si] <= 0) next
    sample_vec <- NULL
    for (fi in which(file_sample_idx == si)) {
      vec <- as.numeric(per_file_rles[[fi]])
      vec[is.na(vec)] <- 0
      sample_vec <- if (is.null(sample_vec)) vec else sample_vec + vec
    }
    sample_vec <- sample_vec / cpm_factors[si]
    acc <- if (length(acc) == 0) sample_vec else acc + sample_vec
  }
  if (length(acc) == 0) next
  cov_vec <- acc / length(samples)
  meanCov <- Rle(cov_vec)
  above <- meanCov >= opt$cutoff
  if (sum(runValue(above)) == 0) next

  # IRanges of TRUE-runs: take the run-encoded ranges, then keep only the
  # ones whose run value is TRUE.
  ir <- ranges(above)[as.logical(runValue(above))]
  if (opt$maxregiongap > 0L) {
    ir <- reduce(ir, min.gapwidth = opt$maxregiongap + 1L)
  }
  if (opt$`min-length` > 1L) {
    ir <- ir[width(ir) >= opt$`min-length`]
  }
  if (length(ir) == 0) next

  scores <- viewMeans(Views(meanCov, ir))
  gr <- GRanges(seqnames = chrom, ranges = ir, strand = "*")
  gr$value <- scores
  regions_per_chrom[[chrom]] <- gr
}

dir.create(dirname(opt$`out-gtf`), recursive = TRUE, showWarnings = FALSE)

if (length(regions_per_chrom) == 0) {
  warning("[run_derfinder] no expressed regions found; writing empty GTF")
  writeLines("# derfinder produced no regions", opt$`out-gtf`)
  quit(status = 0)
}
all_regs <- unlist(GRangesList(regions_per_chrom))
all_regs$gene_id <- paste0("gene", seq_along(all_regs))
all_regs$transcript_id <- paste0("tx", seq_along(all_regs))

# Emit a minimal GTF with one gene + one transcript + one exon per region.
con <- file(opt$`out-gtf`, "w")
on.exit(close(con))
cat("# derfinder regions\n", file = con)
cat(sprintf("# cutoff=%g CPM, min_length=%d bp, maxregiongap=%d bp, chromosomes=%s\n",
            opt$cutoff, opt$`min-length`, opt$maxregiongap,
            paste(chroms, collapse = ",")), file = con)
for (i in seq_along(all_regs)) {
  reg <- all_regs[i]
  chrom <- as.character(seqnames(reg))
  start <- start(reg)
  end <- end(reg)
  strand_chr <- as.character(strand(reg))
  if (strand_chr == "*") strand_chr <- "."
  attrs_gene <- sprintf('gene_id "%s"; gene_name "%s_dergene";',
                        reg$gene_id, reg$gene_id)
  attrs_tx   <- sprintf('gene_id "%s"; transcript_id "%s";',
                        reg$gene_id, reg$transcript_id)
  attrs_exon <- sprintf('gene_id "%s"; transcript_id "%s"; exon_number "1";',
                        reg$gene_id, reg$transcript_id)
  score <- sprintf("%.4f", reg$value[1])
  cat(sprintf("%s\tderfinder\tgene\t%d\t%d\t%s\t%s\t.\t%s\n",
              chrom, start, end, score, strand_chr, attrs_gene), file = con)
  cat(sprintf("%s\tderfinder\ttranscript\t%d\t%d\t%s\t%s\t.\t%s\n",
              chrom, start, end, score, strand_chr, attrs_tx), file = con)
  cat(sprintf("%s\tderfinder\texon\t%d\t%d\t%s\t%s\t.\t%s\n",
              chrom, start, end, score, strand_chr, attrs_exon), file = con)
}
message("[run_derfinder] wrote ", length(all_regs), " regions to ", opt$`out-gtf`)
