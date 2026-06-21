#!/usr/bin/env Rscript
# GTEx structural-concordance heatmap for Main Figure 2 panel G, rebuilt as a
# standalone vector/raster figure from the per-sub-group fastder catalogs. Mirrors
# the Figure 1 chunk of gtex_concordance.Rmd, but drops the per-sub-group row and
# column names: the tissue colour bars and the legend carry the grouping, so the
# 32 sample labels are not needed. Pairwise similarity is the base-pair Jaccard of
# the ER exon segmentations; sub-groups are clustered on 1 - Jaccard.
#
# Usage: make_concordance_heatmap.R [out_basename]
#   out_basename defaults to fig_gtex_concordance (writes .png and .pdf in FIG_DIR).

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(ComplexHeatmap)
  library(circlize)
})

args <- commandArgs(trailingOnly = TRUE)
out_base <- if (length(args) >= 1) args[[1]] else "fig_gtex_concordance"

RESULTS_ROOT <- Sys.getenv("FASTDER_RESULTS_ROOT",
                           "/home/imallona/src/writing_fastder/barbara_results/results")
FIG_DIR <- Sys.getenv("FASTDER_FIG_DIR", "/home/imallona/src/writing_fastder/figures")
CONFIG <- "config_gtex_concordance"

tissue_levels <- c("brain", "heart", "muscle", "blood")
# Tissue palette identical to panel I (panel_novel) so the figure reads as one.
tissue_cols <- c(brain = "#6a51a3", heart = "#cb181d",
                 muscle = "#41ab5d", blood = "#2171b5")

gtf_paths <- Sys.glob(file.path(RESULTS_ROOT, CONFIG, "fastder", "*",
                                "reference", "mc1.0", "gffcompare.annotated.gtf"))
if (length(gtf_paths) == 0)
  stop("no sub-group GTFs under ", file.path(RESULTS_ROOT, CONFIG, "fastder"))

subgroups <- basename(dirname(dirname(dirname(gtf_paths))))
tissue_of <- sub("_[0-9]+$", "", subgroups)
keep <- tissue_of %in% tissue_levels
gtf_paths <- gtf_paths[keep]; subgroups <- subgroups[keep]; tissue_of <- tissue_of[keep]
# Tissue-then-subgroup order keeps the colour bars in blocks before clustering.
ord <- order(factor(tissue_of, levels = tissue_levels), subgroups)
gtf_paths <- gtf_paths[ord]; subgroups <- subgroups[ord]; tissue_of <- tissue_of[ord]

load_exons <- function(path) {
  gr <- tryCatch(import(path), error = function(e) GRanges())
  if (length(gr) == 0) return(GRanges())
  if ("type" %in% colnames(mcols(gr))) {
    ex <- gr[!is.na(gr$type) & gr$type == "exon"]
    if (length(ex) > 0) gr <- ex
  }
  reduce(granges(gr), ignore.strand = TRUE)
}
exons <- lapply(gtf_paths, load_exons)

jaccard_bp <- function(a, b) {
  uni <- sum(as.numeric(width(union(a, b, ignore.strand = TRUE))))
  if (uni == 0) return(NA_real_)
  sum(as.numeric(width(intersect(a, b, ignore.strand = TRUE)))) / uni
}
n <- length(exons)
sim <- matrix(NA_real_, n, n, dimnames = list(subgroups, subgroups))
for (i in seq_len(n)) for (j in i:n) {
  s <- jaccard_bp(exons[[i]], exons[[j]])
  sim[i, j] <- s; sim[j, i] <- s
}

ann_df <- data.frame(tissue = factor(tissue_of, levels = tissue_levels))
top_ann <- HeatmapAnnotation(df = ann_df, col = list(tissue = tissue_cols),
                             show_annotation_name = FALSE)
left_ann <- rowAnnotation(df = ann_df, col = list(tissue = tissue_cols),
                          show_annotation_name = FALSE, show_legend = FALSE)
clust <- hclust(as.dist(1 - sim))

ht <- Heatmap(sim, name = "Jaccard",
              col = colorRamp2(c(0, 0.5, 1), c("white", "#9ecae1", "#08306b")),
              cluster_rows = clust, cluster_columns = clust,
              top_annotation = top_ann, left_annotation = left_ann,
              show_row_names = FALSE, show_column_names = FALSE,
              show_row_dend = FALSE,
              heatmap_legend_param = list(labels_gp = gpar(fontsize = 9),
                                          title_gp = gpar(fontsize = 9)))

png_out <- file.path(FIG_DIR, paste0(out_base, ".png"))
pdf_out <- file.path(FIG_DIR, paste0(out_base, ".pdf"))
png(png_out, width = 5.2, height = 4.2, units = "in", res = 300)
draw(ht, merge_legend = TRUE)
invisible(dev.off())
cairo_pdf(pdf_out, width = 5.2, height = 4.2)
draw(ht, merge_legend = TRUE)
invisible(dev.off())
cat("wrote", png_out, "and", pdf_out, "\n")
