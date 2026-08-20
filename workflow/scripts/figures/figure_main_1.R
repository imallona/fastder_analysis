# Main Figure 1: tool and benchmark.
# A pipeline schematic, B simulation design, C exon accuracy vs depth,
# D boundary precision vs depth, E boundary-distance CDF, F runtime vs memory.
# Data panels are reproduced verbatim from the report code, restyled.
#
# The multi-exon and strand-concordance panels of the first submission are gone.
# Only fastder emits multi-exon structure or assigns a strand, so the bars put
# every other tool at zero by construction, which is a capability statement
# drawn as a measurement. The manuscript states it as a table instead, and the
# fastder quantities stay in the summary report, which is why panel_multiexon()
# and panel_strand() remain in helpers.R.
#
# groHMM is dropped from C and D for the reason argued at EXON_LEVEL_TOOLS in
# helpers.R. It stays in E.

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[[1]] else "figure_main_1.pdf"
fig_dir <- Sys.getenv("FASTDER_FIG_DIR", "/home/imallona/src/writing_fastder/figures")
bench_dir <- Sys.getenv("FASTDER_BENCH_DIR", "/home/imallona/src/writing_fastder/barbara_results/benchmarks/config_full_simulation")

source(file.path(dirname(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE))), "helpers.R"))

p_pipeline  <- wrap_pdf(file.path(fig_dir, "Figure_1.pdf"), width_fill = TRUE)
p_sim       <- wrap_pdf(file.path(fig_dir, "fig_sim_schematic.pdf"))
p_depth     <- panel_depth(which_levels = "Exon", tools = EXON_LEVEL_TOOLS)  # transcript level flat; full version supplementary
p_boundary  <- panel_boundary(tools = EXON_LEVEL_TOOLS)
p_cdf       <- panel_boundary_cdf()
p_speed     <- panel_speed(bench_dir)

# Two full-width schematic rows (A pipeline, B simulation design), then the
# benchmark panels in narrative order: C exon accuracy, D boundary precision,
# E boundary-distance CDF, F speed, each full width now that the two half-width
# capability panels are gone. Add order sets the A-F tags, so it follows the
# order they are cited in the text.
design <- "
AAAAAAAAAAAA
BBBBBBBBBBBB
CCCCCCCCCCCC
DDDDDDDDDDDD
EEEEEEEEEEEE
FFFFFFFFFFFF
"
# Taller rows A (schematic, was cropped) and E (CDF, was crowded).
panels <- p_pipeline + p_sim + p_depth + p_boundary + p_cdf + p_speed +
  plot_layout(design = design, heights = c(4.3, 2.3, 3.4, 2.6, 2.6, 2.4)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

# One dedicated legend (tools as line + symbol, plus strand fill) below the
# panels, so the legend line types match the plotted lines exactly.
fig <- cowplot::plot_grid(panels, make_combined_legend(), ncol = 1, rel_heights = c(1, 0.07))

ggsave(out, fig, width = 8.27, height = 16.0, limitsize = FALSE)
ggsave(sub("\\.pdf$", ".svg", out), fig, width = 8.27, height = 16.0, limitsize = FALSE)
cat("wrote", out, "and", sub("\\.pdf$", ".svg", out), "\n")
