#!/usr/bin/env Rscript
# Exonic Jaccard between the TDP-43 knockdown and control catalogs, the same
# shape metric the GTEx concordance uses (shared exon base pairs over the union),
# at each CPM threshold given. One value per threshold; the catalogs are per
# group, so there is no per-sample clustering to compute.
# Usage: extract_tdp43_jaccard.R <out.csv> <label=manifest.csv> [<label=manifest.csv> ...]

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
out <- args[[1]]
pairs <- args[-1]

load_exons <- function(path) {
  if (!nzchar(path) || !file.exists(path)) return(GRanges())
  gr <- tryCatch(import(path), error = function(e) GRanges())
  if (length(gr) == 0) return(GRanges())
  if ("type" %in% colnames(mcols(gr))) {
    ex <- gr[!is.na(gr$type) & gr$type == "exon"]
    if (length(ex) > 0) gr <- ex
  }
  reduce(granges(gr), ignore.strand = TRUE)
}

jaccard_bp <- function(a, b) {
  inter <- sum(as.numeric(width(GenomicRanges::intersect(a, b, ignore.strand = TRUE))))
  uni <- sum(as.numeric(width(GenomicRanges::union(a, b, ignore.strand = TRUE))))
  if (uni == 0) 0 else inter / uni
}

rows <- lapply(pairs, function(p) {
  label <- sub("=.*$", "", p)
  manifest_csv <- sub("^[^=]*=", "", p)
  manifest <- read.csv(manifest_csv, stringsAsFactors = FALSE)
  kd <- load_exons(unique(manifest$fastder_gtf[manifest$group == "knockdown"])[1])
  ct <- load_exons(unique(manifest$fastder_gtf[manifest$group == "control"])[1])
  data.frame(cpm = label, jaccard = jaccard_bp(kd, ct))
})

write.csv(do.call(rbind, rows), out, row.names = FALSE, quote = FALSE)
