#!/usr/bin/env Rscript
# Novel expressed-region exons per tissue: the union of a tissue's sub-group ER
# exons that overlap no reference exon, matching the gtex_concordance report.
# REF_GTF must be genome-wide, so build figures under config_gtex_concordance.
# Usage: extract_novel_exons.R <fastder_dir> <reference_gtf> <out.csv>

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = TRUE)
fastder_dir <- args[[1]]
reference_gtf <- args[[2]]
out <- args[[3]]

# Disjoint ER exon intervals, strand ignored (recount3 coverage is unstranded).
load_exons <- function(path) {
  if (!file.exists(path)) return(GRanges())
  gr <- tryCatch(import(path), error = function(e) GRanges())
  if (length(gr) == 0) return(GRanges())
  if ("type" %in% colnames(mcols(gr))) {
    ex <- gr[!is.na(gr$type) & gr$type == "exon"]
    if (length(ex) > 0) gr <- ex
  }
  reduce(granges(gr), ignore.strand = TRUE)
}

gtfs <- Sys.glob(file.path(fastder_dir, "*", "reference", "mc1.0",
                           "gffcompare.annotated.gtf"))
subgroups <- basename(dirname(dirname(dirname(gtfs))))  # <tissue>_<n>
tissue_of <- sub("_[0-9]+$", "", subgroups)
tissue_levels <- intersect(c("brain", "heart", "muscle", "blood"), unique(tissue_of))

exons_by_sub <- lapply(gtfs, load_exons)

ann <- import(reference_gtf)
ann <- ann[!is.na(ann$type) & ann$type == "exon"]
seqlevels(ann) <- paste0("chr", sub("^chr", "", seqlevels(ann)))
ref_exons <- reduce(granges(ann), ignore.strand = TRUE)

rows <- lapply(tissue_levels, function(t) {
  parts <- exons_by_sub[tissue_of == t]
  parts <- parts[vapply(parts, length, integer(1)) > 0]
  ex <- if (length(parts))
    reduce(do.call(c, unname(parts)), ignore.strand = TRUE) else GRanges()
  total <- length(ex)
  novel <- if (total) sum(!overlapsAny(ex, ref_exons, ignore.strand = TRUE)) else 0L
  data.frame(tissue = t, novel = novel,
             pct = if (total) round(100 * novel / total) else 0L)
})

write.csv(do.call(rbind, rows), out, row.names = FALSE, quote = FALSE)
