# Main Figure 2: fastder on recount3 data, two worked examples in order.
#   TDP-43: A design, B STMN2 cryptic-exon track (large), C sample similarity,
#           D-F call structure (exons per region, exonic length, ER score)
#   GTEx:   G design, H concordance heatmap (large), I troponin marker loci,
#           J novel ER exons, K exon precision, L transcript precision,
#           M exons per region
# Harmonised cell sizes, sensible per-panel aspect ratios; the user fine-tunes
# placement afterwards. Quantitative panels are vector; track and heatmaps are
# wrapped grid objects.

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[[1]] else "figure_main_2.pdf"
fig_dir <- Sys.getenv("FASTDER_FIG_DIR", "/home/imallona/src/writing_fastder/figures")
TDP <- "config_klim_2019_tdp43_recount3"
GTX <- "config_gtex_comparison"

source(file.path(dirname(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE))), "helpers.R"))

wrap_png <- function(name) {
  wrap_elements(full = grid::rasterGrob(
    png::readPNG(file.path(fig_dir, name)), interpolate = TRUE))
}

A <- wrap_pdf(file.path(fig_dir, "fig_tdp43_scheme.pdf"))
B <- wrap_png("fig_tdp43_stmn2.png")
C <- wrap_png("fig_tdp43_similarity.png")
D <- panel_gtexcmp_exons(TDP)
E <- panel_gtexcmp_length(TDP)
G <- wrap_pdf(file.path(fig_dir, "fig_gtex_scheme.pdf"))
H <- wrap_png("fig_gtex_concordance.png")
I <- panel_marker_loci()
J <- panel_novel()
K <- panel_gtexcmp_precision(GTX, level = "exon")
L <- panel_gtexcmp_precision(GTX, level = "transcript")
M <- panel_gtexcmp_exons(GTX)

design <- "
AAABBBBBBBBB
CCCCCCDDDEEE
FFFGGGGGGGGG
HHHHHHHHHHHH
IIIJJJKKKLLL
"
fig <- A + B + C + D + E + G + H + I + J + K + L + M +
  plot_layout(design = design, heights = c(3.8, 2.0, 2.6, 1.2, 1.8)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))

ggsave(out, fig, width = 8.27, height = 12.8, limitsize = FALSE)
ggsave(sub("\\.pdf$", ".svg", out), fig, width = 8.27, height = 12.8, limitsize = FALSE)
cat("wrote", out, "and", sub("\\.pdf$", ".svg", out), "\n")
