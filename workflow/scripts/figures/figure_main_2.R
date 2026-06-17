# Main Figure 2: fastder on recount3 data, two worked examples in order.
#   TDP-43: design, STMN2 cryptic-exon track, call-structure detail.
#   GTEx:   design, concordance heatmap, marker loci, novel exons, precision.
# Tags follow panel order. Quantitative panels are vector; track and heatmaps
# are wrapped grid objects.

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

# Variable names skip F so they never shadow R's FALSE; the A-L panel tags come
# from the order panels are added, not from these names.
A <- wrap_pdf(file.path(fig_dir, "fig_tdp43_scheme.pdf"))  # TDP-43 design
B <- wrap_png("fig_tdp43_stmn2.png")                       # STMN2 track
D <- panel_gtexcmp_exons(TDP)                              # exons/region, TDP run
E <- panel_gtexcmp_length(TDP)                             # exonic length, TDP run
G <- wrap_pdf(file.path(fig_dir, "fig_gtex_scheme.pdf"))   # GTEx design
H <- wrap_png("fig_gtex_concordance.png")                  # concordance heatmap
I <- panel_marker_loci()                                   # troponin loci
J <- panel_novel()                                         # novel ER exons
K <- panel_gtexcmp_precision(GTX, level = "exon")
L <- panel_gtexcmp_precision(GTX, level = "transcript")
M <- panel_gtexcmp_exons(GTX)                              # exons/region, chr19 run

# Areas are alphabetical by panel order: A scheme, B STMN2, C exons-TDP,
# D length-TDP, E GTEx design, F concordance, G marker loci, H novel, I-K precision/exons.
design <- "
AAABBBBBBBBB
CCCCCCDDDDDD
EEEFFFFFFHHH
GGGGGGGGGGGG
IIIIJJJJKKKK
"
fig <- A + B + D + E + G + H + I + J + K + L + M +
  plot_layout(design = design, heights = c(3.6, 2.0, 3.6, 1.6, 1.9)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))

ggsave(out, fig, width = 8.27, height = 13.2, limitsize = FALSE)
ggsave(sub("\\.pdf$", ".svg", out), fig, width = 8.27, height = 13.2, limitsize = FALSE)
cat("wrote", out, "and", sub("\\.pdf$", ".svg", out), "\n")
