# Main Figure 2: two recount3 worked examples.
#   TDP-43: design, STMN2 track, novel ERs per group.
#   GTEx:   design, concordance heatmap, marker loci, novel exons, precision.

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

A  <- wrap_pdf(file.path(fig_dir, "fig_tdp43_scheme.pdf"))
B  <- wrap_pdf(file.path(fig_dir, "fig_tdp43_stmn2.pdf"))   # vector track, replaces the raster PNG
# C and D are negative results: no genome-wide knockdown/control separation.
C  <- panel_novel_tdp43()
Dd <- panel_tdp43_boundary_dist()
Ej <- panel_tdp43_jaccard()
Fg <- wrap_pdf(file.path(fig_dir, "fig_gtex_scheme.pdf"))
Gc <- wrap_png("fig_gtex_concordance.png")
Hm <- panel_marker_loci()
In <- panel_novel()
Jp <- panel_gtexcmp_precision(GTX, level = "exon")
Kp <- panel_gtexcmp_precision(GTX, level = "transcript")

# Patchwork assigns panels to areas alphabetically. The three small TDP panels
# (C, D, E) sit in one row. A4 width; nudge in the SVG.
design <- "
AAAABBBBBBBB
CCCCDDDDEEEE
FFFFGGGGGGGG
HHHHHIIIIIII
HHHHHJJJJJJJ
HHHHHKKKKKKK
"
fig <- A + B + C + Dd + Ej + Fg + Gc + Hm + In + Jp + Kp +
  plot_layout(design = design, heights = c(2.2, 1.6, 3.0, 1.3, 1.3, 1.3)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11))

ggsave(out, fig, width = 8.27, height = 11, limitsize = FALSE)
ggsave(sub("\\.pdf$", ".svg", out), fig, width = 8.27, height = 11, limitsize = FALSE)
cat("wrote", out, "and", sub("\\.pdf$", ".svg", out), "\n")
