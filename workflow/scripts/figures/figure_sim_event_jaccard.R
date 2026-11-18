# Median best exonic Jaccard of fastder per event class: range over the
# parameter grid per min_coverage, and the reference combination in black.

args <- commandArgs(trailingOnly = TRUE)
fig_dir <- Sys.getenv("FASTDER_FIG_DIR", file.path("results", "figures"))
out <- if (length(args) >= 1) args[[1]] else file.path(fig_dir, "fig_sim_event_jaccard.pdf")

source(file.path(dirname(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE))), "helpers.R"))
BASE_SIZE <- 8

reference_id <- sprintf("mc%s_ml10_pt5_ns0", format(REFERENCE_MIN_COVERAGE, trim = TRUE))
medians <- read_result("config_full_simulation", "fuzzy_jaccard.csv") %>%
  filter(tool == "fastder") %>%
  default_grid() %>%
  mutate(jaccard = as.numeric(jaccard),
         min_coverage = paste(extract_mc(param_id), "CPM")) %>%
  group_by(param_id, min_coverage, scenario, sample) %>%
  summarise(median_jaccard = median(jaccard), .groups = "drop") %>%
  mutate(scenario = relabel_scenario(scenario), event = relabel_event(sample),
         reference = param_id == reference_id)
save_panel_data(medians, "panel_event_jaccard")

ranges <- medians %>%
  group_by(min_coverage, scenario, event) %>%
  summarise(low = min(median_jaccard), high = max(median_jaccard), .groups = "drop")
dodge <- position_dodge(width = 0.5)

p <- ggplot(ranges, aes(event, colour = min_coverage)) +
  geom_linerange(aes(ymin = low, ymax = high), linewidth = 1.6, position = dodge) +
  geom_point(data = filter(medians, reference), aes(y = median_jaccard),
             colour = "black", size = 1.2, position = position_nudge(x = -0.125)) +
  scale_colour_manual(values = c("0.005 CPM" = "#6baed6", "0.05 CPM" = "#fd8d3c"), name = NULL) +
  facet_wrap(~ scenario, ncol = 1) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Event class", y = "Median best exonic Jaccard") +
  theme_pub()

ggsave(out, p, width = 7, height = 4.6)
cat("wrote", out, "\n")
