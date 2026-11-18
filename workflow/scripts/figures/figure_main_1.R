# Main Figure 1: tool and benchmark, two columns at 160 mm width.
# A pipeline overview, B simulation design, C exon accuracy vs depth,
# D boundary precision vs depth, E boundary-distance CDF, F runtime vs memory.
# groHMM leaves C and D; see EXON_LEVEL_TOOLS in helpers.R.

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1) args[[1]] else "figure_main_1.pdf"
fig_dir <- Sys.getenv("FASTDER_FIG_DIR", file.path("results", "figures"))
bench_dir <- Sys.getenv("FASTDER_BENCH_DIR",
                        file.path("logs", "benchmarks", "config_full_simulation"))
script_dir <- dirname(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE)))

source(file.path(script_dir, "helpers.R"))
BASE_SIZE <- 8

p_pipeline <- wrap_pdf(file.path(fig_dir, "fig_pipeline_schematic.pdf"))
p_sim      <- wrap_pdf(file.path(fig_dir, "fig_sim_schematic.pdf"))
p_depth    <- free(panel_depth(which_levels = "Exon", tools = EXON_LEVEL_TOOLS))
p_boundary <- panel_boundary(tools = EXON_LEVEL_TOOLS)
p_cdf      <- panel_boundary_cdf()
p_speed    <- panel_speed(bench_dir)

design <- "
AAAAAAAAAAAA
BBBBBBBFFFFF
CCCCCCCDDDDD
CCCCCCCEEEEE
"
panels <- p_pipeline + p_sim + p_depth + p_boundary + p_cdf + p_speed +
  plot_layout(design = design, heights = c(1.3, 2.1, 1.6, 1.6)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 10))

fig <- cowplot::plot_grid(panels, make_combined_legend(strand = FALSE),
                          ncol = 1, rel_heights = c(1, 0.04))

width_in <- 160 / 25.4
height_in <- 160 / 25.4
ggsave(out, fig, width = width_in, height = height_in)
ggsave(sub("\\.pdf$", ".svg", out), fig, width = width_in, height = height_in)
cat("wrote", out, "and", sub("\\.pdf$", ".svg", out), "\n")
