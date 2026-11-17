# Median best exonic Jaccard along the tolerance axis of fastder
# (position_tolerance) and derfinder (maxRegionGap), other parameters at the
# reference point.

args <- commandArgs(trailingOnly = TRUE)
fig_dir <- Sys.getenv("FASTDER_FIG_DIR", file.path("results", "figures"))
out <- if (length(args) >= 1) args[[1]] else file.path(fig_dir, "fig_sim_tolerance.pdf")

source(file.path(dirname(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE))), "helpers.R"))
BASE_SIZE <- 8

axis_labels <- c(fastder = "fastder, position_tolerance", derfinder = "derfinder, maxRegionGap")
sweep <- read_result("config_full_simulation", "fuzzy_jaccard.csv") %>%
  filter(tool %in% names(axis_labels)) %>%
  default_grid() %>% at_reference_threshold() %>%
  filter(tool != "fastder" | grepl("_ml10_", param_id)) %>%
  mutate(jaccard = as.numeric(jaccard), tolerance = extract_pt(param_id)) %>%
  group_by(tool, scenario, tolerance) %>%
  summarise(median_jaccard = median(jaccard), .groups = "drop") %>%
  mutate(scenario = relabel_scenario(scenario),
         axis = factor(axis_labels[tool], levels = axis_labels))
save_panel_data(sweep, "panel_tolerance")

p <- ggplot(sweep, aes(tolerance, median_jaccard, colour = tool, shape = tool)) +
  geom_line(linewidth = 0.6) + geom_point(size = 2) +
  scale_colour_manual(values = tool_palette, guide = "none") +
  scale_shape_manual(values = tool_shapes, guide = "none") +
  scale_x_continuous(breaks = c(0, 5, 10, 20)) +
  facet_grid(scenario ~ axis, labeller = labeller(scenario = label_wrap_gen(12))) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Tolerance (bp)", y = "Median best exonic Jaccard") +
  theme_pub_square()

ggsave(out, p, width = 4.6, height = 4.4)
cat("wrote", out, "\n")
