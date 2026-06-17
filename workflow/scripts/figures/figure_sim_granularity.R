# Predicted-transcript granularity with the simulated truth marked as a cross.
# Reads config_full_simulation; rebuilds from the local mirror.

args <- commandArgs(trailingOnly = TRUE)
fig_dir <- Sys.getenv("FASTDER_FIG_DIR", "/home/imallona/src/writing_fastder/figures")
out <- if (length(args) >= 1) args[[1]] else file.path(fig_dir, "fig_sim_granularity.pdf")

source(file.path(dirname(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE))), "helpers.R"))
suppressPackageStartupMessages(library(ggrepel))

config <- "config_full_simulation"
jac <- read_result(config, "fuzzy_jaccard.csv")
pids <- best_pids(jac)
fj <- pick_tool_param(jac, pids) %>%
  mutate(fastder_bp = as.numeric(fastder_bp), scenario = relabel_scenario(scenario))
granularity <- fj %>%
  filter(!is.na(fastder_transcript), fastder_transcript != "") %>%
  distinct(tool, scenario, fastder_transcript, fastder_bp) %>%
  group_by(tool, scenario) %>%
  summarise(n_transcripts = n(), median_bp = median(fastder_bp, na.rm = TRUE), .groups = "drop")
truth <- read_csv(file.path(RESULTS_ROOT, config, "truth_chain_stats.csv"), show_col_types = FALSE) %>%
  mutate(scenario = relabel_scenario(scenario)) %>%
  group_by(scenario) %>%
  summarise(n_transcripts = n(), median_bp = median(total_exon_length, na.rm = TRUE), .groups = "drop")

p <- ggplot(granularity, aes(n_transcripts, median_bp, colour = tool, shape = tool)) +
  geom_point(size = 4, alpha = 0.85) +
  geom_text_repel(aes(label = tool), size = 3, colour = "grey20", box.padding = 0.5, max.overlaps = Inf) +
  scale_colour_manual(values = tool_palette, name = NULL) +
  scale_shape_manual(values = tool_shapes, name = NULL) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~ scenario, ncol = 2) +
  labs(x = "Number of predicted transcripts (log)", y = "Median predicted transcript length, bp (log)") +
  theme_pub() + theme(legend.position = "none") +
  geom_point(data = truth, aes(n_transcripts, median_bp), inherit.aes = FALSE,
             shape = 4, size = 3, stroke = 0.9, colour = "black") +
  geom_text(data = truth, aes(n_transcripts, median_bp, label = "truth"), inherit.aes = FALSE,
            vjust = -1.1, size = 2.6, colour = "black")

ggsave(out, p, width = 7, height = 5)
ggsave(sub("\\.pdf$", ".png", out), p, width = 7, height = 5, dpi = 300)
ggsave(sub("\\.pdf$", ".svg", out), p, width = 7, height = 5)
cat("wrote", out, "(+ png, svg)\n")
