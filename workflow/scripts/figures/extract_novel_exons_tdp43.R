#!/usr/bin/env Rscript
# Novel ER exons per group (knockdown vs control) for the TDP-43 run: fastder ER
# exons overlapping no Ensembl exon. Novel definition as in extract_novel_exons.R,
# but the unit is the group and the GTF comes from the manifest fastder_gtf
# column. Reference must cover the called chromosomes.
# Usage: extract_novel_exons_tdp43.R <recount3_manifest.csv> <reference_gtf> <out.csv>

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
manifest_csv <- args[[1]]
reference_gtf <- args[[2]]
out <- args[[3]]

# Disjoint ER exon intervals, strand ignored (recount3 coverage is unstranded).
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

manifest <- read.csv(manifest_csv, stringsAsFactors = FALSE)

ann <- import(reference_gtf)
ann <- ann[!is.na(ann$type) & ann$type == "exon"]
seqlevels(ann) <- paste0("chr", sub("^chr", "", seqlevels(ann)))
ref_exons <- reduce(granges(ann), ignore.strand = TRUE)

# Knockdown first so the panel reads kd then control.
groups <- intersect(c("knockdown", "control"), unique(manifest$group))
rows <- lapply(groups, function(g) {
  gtf <- unique(manifest$fastder_gtf[manifest$group == g])
  ex <- load_exons(if (length(gtf)) gtf[1] else "")
  total <- length(ex)
  novel <- if (total) sum(!overlapsAny(ex, ref_exons, ignore.strand = TRUE)) else 0L
  data.frame(group = g, novel = novel,
             pct = if (total) round(100 * novel / total) else 0L)
})

write.csv(do.call(rbind, rows), out, row.names = FALSE, quote = FALSE)
