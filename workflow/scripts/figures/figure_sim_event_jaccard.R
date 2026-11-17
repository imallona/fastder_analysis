# Median best exonic Jaccard of fastder per event class, one line per
# parameter combination of the grid; the reference combination is drawn black.

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

p <- ggplot(medians, aes(event, median_jaccard, group = param_id)) +
  geom_line(data = filter(medians, !reference), aes(colour = min_coverage),
            linewidth = 0.35, alpha = 0.7) +
  geom_line(data = filter(medians, reference), colour = "black", linewidth = 0.8) +
  geom_point(data = filter(medians, reference), colour = "black", size = 1.4) +
  scale_colour_manual(values = c("0.005 CPM" = "#9ecae1", "0.05 CPM" = "#fdae6b"), name = NULL) +
  facet_wrap(~ scenario, ncol = 1) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Event class", y = "Median best exonic Jaccard") +
  theme_pub()

ggsave(out, p, width = 7, height = 4.6)
cat("wrote", out, "\n")
