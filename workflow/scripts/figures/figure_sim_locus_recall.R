# Locus recall at 50% exonic overlap per tool and event class, at the
# parameters of the reported numbers.

args <- commandArgs(trailingOnly = TRUE)
fig_dir <- Sys.getenv("FASTDER_FIG_DIR", file.path("results", "figures"))
out <- if (length(args) >= 1) args[[1]] else file.path(fig_dir, "fig_sim_locus_recall.pdf")

source(file.path(dirname(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE))), "helpers.R"))
BASE_SIZE <- 8

config <- "config_full_simulation"
pids <- reference_pids(read_result(config, "fuzzy_jaccard.csv"))
recall <- read_result(config, "fuzzy_locus_recall.csv") %>%
  pick_tool_param(pids) %>%
  mutate(threshold = as.numeric(threshold), recall = 100 * as.numeric(recall)) %>%
  filter(threshold == 0.5) %>%
  mutate(scenario = relabel_scenario(scenario), event = relabel_event(sample))
save_panel_data(recall, "panel_locus_recall")

p <- ggplot(recall, aes(event, recall, colour = tool, shape = tool)) +
  geom_point(size = 2, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = tool_palette, labels = tool_labels, name = NULL) +
  scale_shape_manual(values = tool_shapes, labels = tool_labels, name = NULL) +
  facet_wrap(~ scenario, ncol = 1) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Event class", y = "Locus recall (%)") +
  theme_pub()

ggsave(out, p, width = 7, height = 4.6)
cat("wrote", out, "\n")
