# Regions per sample that match a simulated transcript, and their median exonic
# length, per tool and for the truth, at the parameters of the reported numbers.

args <- commandArgs(trailingOnly = TRUE)
fig_dir <- Sys.getenv("FASTDER_FIG_DIR", file.path("results", "figures"))
out <- if (length(args) >= 1) args[[1]] else file.path(fig_dir, "fig_sim_granularity.pdf")

source(file.path(dirname(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE))), "helpers.R"))
BASE_SIZE <- 8

config <- "config_full_simulation"
jac <- read_result(config, "fuzzy_jaccard.csv")
fj <- pick_tool_param(jac, reference_pids(jac)) %>%
  mutate(fastder_bp = as.numeric(fastder_bp), scenario = relabel_scenario(scenario))
granularity <- fj %>%
  filter(!is.na(fastder_transcript), fastder_transcript != "") %>%
  distinct(tool, scenario, sample, fastder_transcript, fastder_bp) %>%
  group_by(tool, scenario) %>%
  summarise(n_regions = n() / n_distinct(sample), median_bp = median(fastder_bp, na.rm = TRUE),
            .groups = "drop")
truth <- read_csv(file.path(RESULTS_ROOT, config, "truth_chain_stats.csv"), show_col_types = FALSE) %>%
  mutate(scenario = relabel_scenario(scenario)) %>%
  group_by(scenario) %>%
  summarise(n_regions = n() / n_distinct(sample), median_bp = median(total_exon_length, na.rm = TRUE),
            .groups = "drop")
save_panel_data(granularity, "panel_granularity")

p <- ggplot(granularity, aes(n_regions, median_bp, colour = tool, shape = tool)) +
  geom_point(size = 2.4) +
  geom_point(data = truth, aes(n_regions, median_bp), inherit.aes = FALSE,
             shape = 4, size = 2.4, stroke = 0.8, colour = "black") +
  scale_colour_manual(values = c(tool_palette, truth = "black"), labels = tool_labels, name = NULL) +
  scale_shape_manual(values = tool_shapes, labels = tool_labels, name = NULL) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~ scenario, ncol = 2) +
  labs(x = "Matched regions per sample (log scale)", y = "Median exonic length, bp (log scale)") +
  theme_pub_square()

ggsave(out, p, width = 5.5, height = 3.2)
cat("wrote", out, "\n")
