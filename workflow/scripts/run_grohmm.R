#!/usr/bin/env Rscript
# HMM-based expressed-region caller, used as a third baseline against
# fastder and the derfinder / megadepth-baseline pair.
#
# Reads every per-sample BigWig in --bigwig-dir, normalises each sample's
# coverage to CPM (counts per million) using the same library_size formula
# fastder and run_derfinder.R apply (Sum width * value over the user's
# chromosomes), then summarises mean per 50 bp window with the kent
# bigWigAverageOverBed utility. Per-window CPMs are averaged across samples
# into one vector per chromosome, integer-scaled, and fed to groHMM's HMM
# via detectTranscripts. groHMM requires both Fp and Fm; since recount3
# coverage is unstranded we feed the same vector as both and deduplicate
# the per-strand calls after unstranding. The output transcripts are
# written as unstranded GTF, matching the derfinder / megadepth-baseline
# outputs we are comparing against.
#
# Uses kent's bigWigAverageOverBed for the per-window summarisation rather
# than rtracklayer::import.bw: the binary avoids R object materialisation and
# is several orders of magnitude faster on whole-chromosome scans.
#
# Parameter equivalence with the other tools:
#   --ltprobb     <-> groHMM LtProbB (log-prob of remaining in the transcribed
#                     state; more negative = stricter calling; default -200).
#   --uts         <-> groHMM UTS (variance of the untranscribed state; default 5).
#   --window-size <-> groHMM's binning width; default 50 bp.
#   --min-length  <-> post-filter on the called intervals (bp).
#   --chromosomes <-> scope used for library_size and for window construction.
suppressPackageStartupMessages({
  library(optparse)
  library(groHMM)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})

opt_list <- list(
  make_option("--bigwig-dir", type = "character", help = "Directory holding per-sample .bw files"),
  make_option("--out-gtf", type = "character", help = "Path of the GTF to write"),
  make_option("--ltprobb", type = "double", default = -200,
              help = "groHMM LtProbB (log-prob of staying in the transcribed state)"),
  make_option("--uts", type = "double", default = 5,
              help = "groHMM UTS (variance of the untranscribed state)"),
  make_option("--window-size", type = "integer", default = 50,
              help = "Window width in bp (groHMM default is 50)"),
  make_option("--min-length", type = "integer", default = 10,
              help = "Drop called intervals shorter than this many bases"),
  make_option("--chromosomes", type = "character", default = NULL,
              help = "Space-separated list of chromosomes to analyse"),
  make_option("--count-scale", type = "integer", default = 100,
              help = "Integer scale applied to mean CPM before HMM fit; controls dynamic range")
)
parser <- OptionParser(option_list = opt_list)
arg_strings <- commandArgs(trailingOnly = TRUE)
# --chromosomes is a space-separated multi-value flag, like in run_derfinder.R.
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

# TODO(deprecate): the .plus.bw / .minus.bw branch supports the workflow's
# stranded=true hack (see config.yaml). No paper figure uses it, but
# removing it needs the original stranded-path contributor's agreement.
# Until then we keep the branch and collapse both strands per sample below.
bw_files <- list.files(opt$`bigwig-dir`, pattern = "\\.(all|plus|minus)\\.bw$",
                       full.names = TRUE)
if (length(bw_files) == 0) {
  stop("[run_grohmm] No .bw files in ", opt$`bigwig-dir`)
}
sample_of <- sub("\\.(all|plus|minus)\\.bw$", "", basename(bw_files))
samples <- split(unname(bw_files), sample_of)
sample_ids <- names(samples)
flat_files <- unlist(samples, use.names = FALSE)
file_sample_idx <- rep(seq_along(samples), lengths(samples))

# Chromosome sizes come from the first BigWig via kent's bigWigInfo. The
# output is one "chrom size" line per chromosome.
read_chrom_sizes <- function(bw_path) {
  raw <- system2("bigWigInfo", c("-chroms", bw_path), stdout = TRUE)
  raw <- raw[grepl("^\\s+", raw)]
  fields <- do.call(rbind, strsplit(trimws(raw), "\\s+"))
  out <- as.integer(fields[, 3])
  names(out) <- fields[, 1]
  out
}

all_sizes <- read_chrom_sizes(flat_files[1])
if (length(chroms) == 0) chroms <- names(all_sizes)
chroms <- intersect(chroms, names(all_sizes))
if (length(chroms) == 0) {
  stop("[run_grohmm] None of the requested chromosomes are in the BigWigs")
}
sizes <- all_sizes[chroms]
message("[run_grohmm] processing ", length(chroms), " chromosomes; ",
        length(samples), " samples")

# 50 bp window BED. Emit each window as a single line so bigWigAverageOverBed
# returns one mean per window. The name column carries chrom and 0-based start
# so we can pivot the TSV back into per-chromosome vectors.
windows_bed <- tempfile(fileext = ".bed")
window_size <- as.integer(opt$`window-size`)
con <- file(windows_bed, "w")
on.exit(close(con), add = TRUE)
for (chrom in chroms) {
  chrom_len <- sizes[[chrom]]
  starts <- seq.int(0L, chrom_len - 1L, by = window_size)
  ends <- pmin(starts + window_size, chrom_len)
  names_col <- paste0(chrom, ":", starts)
  cat(sprintf("%s\t%d\t%d\t%s\n", chrom, starts, ends, names_col),
      file = con, sep = "")
}
close(con); on.exit()

# Per-sample, per-chromosome library size (Sum width * value) restricted to
# the requested chromosomes, computed once per file with bigWigInfo. fastder
# and run_derfinder.R both use this exact formula.
compute_library_size_kent <- function(bw_path, target_chroms) {
  raw <- system2("bigWigInfo", c("-chroms", bw_path), stdout = TRUE)
  raw <- raw[grepl("^\\s+", raw)]
  fields <- do.call(rbind, strsplit(trimws(raw), "\\s+"))
  available_chroms <- fields[, 1]
  keep <- intersect(target_chroms, available_chroms)
  if (length(keep) == 0) return(0)
  # bigWigAverageOverBed with a single per-chrom interval returns sum = sum of
  # values inside the interval, which equals Sum width * value for that
  # chromosome (after multiplication by interval width = 1 in bedGraph units).
  # We use the per-chrom mean times chrom length to get the same number
  # without writing a chrom-spanning BED.
  total <- 0
  for (chrom in keep) {
    chrom_len <- as.integer(fields[match(chrom, available_chroms), 3])
    bed <- tempfile(fileext = ".bed")
    writeLines(sprintf("%s\t0\t%d\t%s", chrom, chrom_len, chrom), bed)
    out <- system2("bigWigAverageOverBed", c(bw_path, bed, "/dev/stdout"),
                   stdout = TRUE)
    file.remove(bed)
    if (length(out) == 0) next
    parts <- strsplit(out, "\t", fixed = TRUE)[[1]]
    # columns: name, size, covered, sum, mean0, mean
    total <- total + as.numeric(parts[4])
  }
  total
}

message("[run_grohmm] computing library sizes")
lib_sizes <- vapply(samples, function(files)
  sum(vapply(files, compute_library_size_kent, numeric(1),
             target_chroms = chroms)),
  FUN.VALUE = numeric(1))
cpm_factors <- lib_sizes / 1e6
for (i in seq_along(samples)) {
  if (cpm_factors[i] <= 0) {
    message("[run_grohmm] WARN: sample ", sample_ids[i],
            " has empty library_size on the target chromosomes; skipping")
  }
}

# Pre-compute the per-chromosome window count and the index range that each
# chromosome occupies in the global per-window vector returned by
# bigWigAverageOverBed (rows come back in BED order, which is the order we
# emitted: chromosome by chromosome).
windows_per_chrom <- vapply(chroms, function(chrom) {
  chrom_len <- sizes[[chrom]]
  length(seq.int(0L, chrom_len - 1L, by = window_size))
}, integer(1))
chrom_offsets <- c(0L, cumsum(windows_per_chrom)[-length(windows_per_chrom)])
names(chrom_offsets) <- chroms
total_windows <- sum(windows_per_chrom)
mean_cpm_sum <- numeric(total_windows)
n_contributing <- 0L

for (si in seq_along(samples)) {
  if (cpm_factors[si] <= 0) next
  per_file_means <- numeric(total_windows)
  for (fi in which(file_sample_idx == si)) {
    bw_path <- flat_files[fi]
    out_tsv <- tempfile(fileext = ".tsv")
    system2("bigWigAverageOverBed", c(bw_path, windows_bed, out_tsv))
    df <- read.table(out_tsv, sep = "\t", header = FALSE,
                     stringsAsFactors = FALSE,
                     colClasses = c("character", "integer", "integer",
                                    "numeric", "numeric", "numeric"))
    file.remove(out_tsv)
    if (nrow(df) != total_windows) {
      stop("[run_grohmm] bigWigAverageOverBed returned ", nrow(df),
           " rows, expected ", total_windows, " for ", bw_path)
    }
    # column 5 is mean0 (mean coverage including uncovered bases as 0); that
    # is the right denominator for a per-window CPM.
    per_file_means <- per_file_means + df[[5]]
  }
  mean_cpm_sum <- mean_cpm_sum + per_file_means / cpm_factors[si]
  n_contributing <- n_contributing + 1L
}
file.remove(windows_bed); on.exit()

if (n_contributing == 0L) {
  warning("[run_grohmm] no samples contributed CPM; writing empty GTF")
  dir.create(dirname(opt$`out-gtf`), recursive = TRUE, showWarnings = FALSE)
  writeLines("# groHMM produced no regions", opt$`out-gtf`)
  quit(status = 0)
}

mean_cpm <- mean_cpm_sum / n_contributing

# Integer-scale the float CPMs so groHMM's HMM has a count-shaped input. The
# LtProbB / UTS parameters absorb the scale, so this constant only sets dynamic
# range. 100 gives roughly 1 unit per 0.01 CPM, which keeps a CPM=0.05 floor
# at 5 units; that is enough granularity for the HMM to separate from zero.
scaled <- as.integer(round(mean_cpm * opt$`count-scale`))
scaled[scaled < 0L] <- 0L

Fp_list <- vector("list", length(chroms))
names(Fp_list) <- chroms
for (chrom in chroms) {
  off <- chrom_offsets[[chrom]]
  n <- windows_per_chrom[[chrom]]
  Fp_list[[chrom]] <- Rle(scaled[(off + 1L):(off + n)])
}

message("[run_grohmm] fitting HMM (LtProbB = ", opt$ltprobb,
        ", UTS = ", opt$uts, ", windows = ", total_windows, ")")
# detectTranscripts requires both Fp and Fm; recount3 coverage is unstranded
# so we feed the same per-window vector as both. The HMM is deterministic on
# identical input, so + and - produce identical intervals; we deduplicate
# them below after unstranding.
res <- detectTranscripts(Fp = Fp_list, Fm = Fp_list, LtProbB = opt$ltprobb,
                         UTS = opt$uts, threshold = 1)
calls <- res$transcripts
if (is.null(calls) || length(calls) == 0L) {
  warning("[run_grohmm] HMM produced no transcripts; writing empty GTF")
  dir.create(dirname(opt$`out-gtf`), recursive = TRUE, showWarnings = FALSE)
  writeLines("# groHMM produced no regions", opt$`out-gtf`)
  quit(status = 0)
}

# detectTranscripts returns transcripts in window units. Convert back to bp
# and unstrand (recount3 coverage is unstranded; the other two tools also
# output unstranded calls on the same input). Because we fed identical Fp
# and Fm above, each region is returned twice; unique() collapses the pair.
strand(calls) <- "*"
calls <- unique(calls)
calls <- calls[as.character(seqnames(calls)) %in% chroms]
if (opt$`min-length` > 1L) {
  calls <- calls[width(calls) >= opt$`min-length`]
}
if (length(calls) == 0L) {
  warning("[run_grohmm] all calls were filtered by --min-length; writing empty GTF")
  dir.create(dirname(opt$`out-gtf`), recursive = TRUE, showWarnings = FALSE)
  writeLines("# groHMM produced no regions", opt$`out-gtf`)
  quit(status = 0)
}

calls$gene_id <- paste0("gene", seq_along(calls))
calls$transcript_id <- paste0("tx", seq_along(calls))
calls$score <- vapply(seq_along(calls), function(i) {
  chrom <- as.character(seqnames(calls)[i])
  s <- start(calls)[i]; e <- end(calls)[i]
  w_start <- (s - 1L) %/% window_size + 1L
  w_end <- (e - 1L) %/% window_size + 1L
  off <- chrom_offsets[[chrom]]
  mean(mean_cpm[(off + w_start):(off + w_end)])
}, numeric(1))

dir.create(dirname(opt$`out-gtf`), recursive = TRUE, showWarnings = FALSE)
con_out <- file(opt$`out-gtf`, "w")
on.exit(close(con_out), add = TRUE)
cat("# groHMM regions\n", file = con_out)
cat(sprintf("# LtProbB=%g, UTS=%g, window_size=%d bp, min_length=%d bp, chromosomes=%s\n",
            opt$ltprobb, opt$uts, window_size, opt$`min-length`,
            paste(chroms, collapse = ",")), file = con_out)
for (i in seq_along(calls)) {
  chrom <- as.character(seqnames(calls)[i])
  s <- start(calls)[i]; e <- end(calls)[i]
  score <- sprintf("%.4f", calls$score[i])
  attrs_gene <- sprintf('gene_id "%s"; gene_name "%s_grohmmgene";',
                        calls$gene_id[i], calls$gene_id[i])
  attrs_tx <- sprintf('gene_id "%s"; transcript_id "%s";',
                      calls$gene_id[i], calls$transcript_id[i])
  attrs_exon <- sprintf('gene_id "%s"; transcript_id "%s"; exon_number "1";',
                        calls$gene_id[i], calls$transcript_id[i])
  cat(sprintf("%s\tgroHMM\tgene\t%d\t%d\t%s\t.\t.\t%s\n",
              chrom, s, e, score, attrs_gene), file = con_out)
  cat(sprintf("%s\tgroHMM\ttranscript\t%d\t%d\t%s\t.\t.\t%s\n",
              chrom, s, e, score, attrs_tx), file = con_out)
  cat(sprintf("%s\tgroHMM\texon\t%d\t%d\t%s\t.\t.\t%s\n",
              chrom, s, e, score, attrs_exon), file = con_out)
}
message("[run_grohmm] wrote ", length(calls), " regions to ", opt$`out-gtf`)
