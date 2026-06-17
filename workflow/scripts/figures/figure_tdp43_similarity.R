# TDP-43 knockdown vs control structural similarity: exonic base-pair Jaccard of
# the fastder ER segmentations (2x2).
# Needs the per-group gffcompare.annotated.gtf pulled from barbara.

args <- commandArgs(trailingOnly = TRUE)
fig_dir <- Sys.getenv("FASTDER_FIG_DIR", "/home/imallona/src/writing_fastder/figures")
out <- if (length(args) >= 1) args[[1]] else file.path(fig_dir, "fig_tdp43_similarity.pdf")

source(file.path(dirname(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE))), "helpers.R"))
suppressPackageStartupMessages({library(rtracklayer); library(GenomicRanges)})

config <- "config_klim_2019_tdp43_recount3"
groups <- c(knockdown = "knockdown", control = "control")
gtf <- file.path(RESULTS_ROOT, config, "fastder", groups, "reference/mc1.0/gffcompare.annotated.gtf")

load_exons <- function(path) {
  gr <- import(path)
  if ("type" %in% colnames(mcols(gr))) gr <- gr[!is.na(gr$type) & gr$type == "exon"]
  reduce(granges(gr), ignore.strand = TRUE)
}
ex <- lapply(gtf, load_exons)

jaccard_bp <- function(a, b) {
  uni <- sum(as.numeric(width(union(a, b, ignore.strand = TRUE))))
  if (uni == 0) NA_real_ else sum(as.numeric(width(intersect(a, b, ignore.strand = TRUE)))) / uni
}
g <- names(groups)
sim <- outer(seq_along(g), seq_along(g), Vectorize(function(i, j) jaccard_bp(ex[[i]], ex[[j]])))
dimnames(sim) <- list(g, g)

long <- expand.grid(s1 = factor(g, levels = g), s2 = factor(g, levels = rev(g)))
long$jaccard <- mapply(function(a, b) sim[a, b], as.character(long$s1), as.character(long$s2))
p <- ggplot(long, aes(s1, s2, fill = jaccard)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.3f", jaccard)), size = 5, colour = "white") +
  scale_fill_gradientn(colours = c("white", "#9ecae1", "#08306b"), values = c(0, 0.5, 1),
                       limits = c(0, 1), name = "Jaccard") +
  coord_fixed() +
  labs(x = NULL, y = NULL, title = "Expressed-region structural similarity") +
  theme_pub() + theme(legend.position = "right", legend.title = element_text())

ggsave(out, p, width = 5, height = 4)
ggsave(sub("\\.pdf$", ".png", out), p, width = 5, height = 4, dpi = 300)
ggsave(sub("\\.pdf$", ".svg", out), p, width = 5, height = 4)
cat("kd vs control Jaccard:", round(sim["knockdown", "control"], 3), "\n")
cat("wrote", out, "(+ png, svg)\n")
