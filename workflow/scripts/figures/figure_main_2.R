# Main Figure 2: two recount3 examples, two columns at 160 mm width.
#   TDP-43: A design, B novel exons per group, C knockdown against control
#           Jaccard, D STMN2 track.
#   GTEx:   E design, F concordance heatmap, G novel exons per tissue,
#           H marker loci, I exon precision per catalog.

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[[1]] else "figure_main_2.pdf"
fig_dir <- Sys.getenv("FASTDER_FIG_DIR", file.path("results", "figures"))
GTX <- "config_gtex_comparison"

source(file.path(dirname(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE))), "helpers.R"))
BASE_SIZE <- 8

tdp_design   <- wrap_pdf(file.path(fig_dir, "fig_tdp43_scheme.pdf"))
tdp_novel    <- panel_novel_tdp43()
tdp_jaccard  <- panel_tdp43_jaccard()
stmn2_track  <- wrap_pdf(file.path(fig_dir, "fig_tdp43_stmn2.pdf"))
gtex_design  <- wrap_pdf(file.path(fig_dir, "fig_gtex_scheme.pdf"))
concordance  <- wrap_png("fig_gtex_concordance.png")
gtex_novel   <- panel_novel()
marker_loci  <- panel_marker_loci()
precision    <- panel_gtexcmp_precision(GTX, level = "exon")

design <- "
AAAABBBBCCCC
DDDDDDDDDDDD
EEEEFFFFFFFF
GGGGFFFFFFFF
HHHHHHHHIIII
"
fig <- tdp_design + tdp_novel + tdp_jaccard + stmn2_track + gtex_design +
  concordance + gtex_novel + marker_loci + precision +
  plot_layout(design = design, heights = c(1.5, 3.0, 1.2, 1.6, 1.6)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 10),
        plot.tag.location = "margin")

width_in <- 160 / 25.4
height_in <- 195 / 25.4
ggsave(out, fig, width = width_in, height = height_in)
ggsave(sub("\\.pdf$", ".svg", out), fig, width = width_in, height = height_in)
cat("wrote", out, "and", sub("\\.pdf$", ".svg", out), "\n")
