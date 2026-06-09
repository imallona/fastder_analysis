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

# Variable names skip F so they never shadow R's FALSE; the A-L panel tags come
# from the order panels are added, not from these names.
A <- wrap_pdf(file.path(fig_dir, "fig_tdp43_scheme.pdf"))  # tag A
B <- wrap_png("fig_tdp43_stmn2.png")                       # tag B  STMN2 track
C <- wrap_png("fig_tdp43_similarity.png")                  # tag C
D <- panel_gtexcmp_exons(TDP)                              # tag D  exons/region, TDP run
E <- panel_gtexcmp_length(TDP)                             # tag E  exonic length, TDP run
G <- wrap_pdf(file.path(fig_dir, "fig_gtex_scheme.pdf"))   # tag F  GTEx design
H <- wrap_png("fig_gtex_concordance.png")                  # tag G  concordance heatmap
I <- panel_marker_loci()                                   # tag H  troponin loci
J <- panel_novel()                                         # tag I  novel ER exons
K <- panel_gtexcmp_precision(GTX, level = "exon")          # tag J
L <- panel_gtexcmp_precision(GTX, level = "transcript")    # tag K
M <- panel_gtexcmp_exons(GTX)                              # tag L  exons/region, chr19 run

# TDP-43 block (rows 1-2): scheme + STMN2 track, then call-structure detail.
# GTEx block (rows 3-5): the square concordance heatmap gets a near-square
# centre block, flanked by the GTEx design (tag F) and the novel-exon bars
# (tag I); then the troponin loci full width; then the chr19 comparison bars.
design <- "
AAABBBBBBBBB
CCCCDDDDEEEE
FFFGGGGGGIII
HHHHHHHHHHHH
JJJJKKKKLLLL
"
fig <- A + B + C + D + E + G + H + I + J + K + L + M +
  plot_layout(design = design, heights = c(3.6, 2.0, 3.6, 1.6, 1.9)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 15))

ggsave(out, fig, width = 8.27, height = 13.2, limitsize = FALSE)
ggsave(sub("\\.pdf$", ".svg", out), fig, width = 8.27, height = 13.2, limitsize = FALSE)
cat("wrote", out, "and", sub("\\.pdf$", ".svg", out), "\n")
