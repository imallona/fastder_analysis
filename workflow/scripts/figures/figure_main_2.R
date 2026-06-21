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

# A and F use native fit placement (no stretch) at identical source size, so the
# two design headers print at the same font size. B is delivered standalone
# (fig_tdp43_stmn2.pdf) and placed into this reserved slot by hand.
A  <- wrap_pdf(file.path(fig_dir, "fig_tdp43_scheme.pdf"))
B  <- panel_placeholder("panel B: STMN2 coverage track\n(place fig_tdp43_stmn2.pdf here)")
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
# (C, D, E) sit in one row. B (STMN2 track) and G (concordance heatmap) each
# span two design rows so they render at double height; A and F sit beside them
# in the top row only, with empty space below. A4 width; nudge in the SVG.
# A (TDP design) and F (GTEx design) get identical 6-col x 2-row cells so the two
# headers match. B's slot is reserved for the standalone track. C/D/E share a row.
design <- "
AAAAAABBBBBB
AAAAAABBBBBB
CCCCDDDDEEEE
FFFFFF######
FFFFFF######
GGGGGGGGGGGG
GGGGGGGGGGGG
HHHHHHHHHHHH
IIIIJJJJKKKK
"
fig <- A + B + C + Dd + Ej + Fg + Gc + Hm + In + Jp + Kp +
  plot_layout(design = design,
              heights = c(1.1, 1.1, 1.7, 1.1, 1.1, 1.9, 1.9, 2.0, 1.9)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 11))

ggsave(out, fig, width = 8.27, height = 13.5, limitsize = FALSE)
ggsave(sub("\\.pdf$", ".svg", out), fig, width = 8.27, height = 13.5, limitsize = FALSE)
cat("wrote", out, "and", sub("\\.pdf$", ".svg", out), "\n")
