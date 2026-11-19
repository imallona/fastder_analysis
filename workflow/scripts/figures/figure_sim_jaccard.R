# Best exonic Jaccard of each simulated transcript against the calls, per tool
# and event class, at the parameters of the reported numbers.

args <- commandArgs(trailingOnly = TRUE)
fig_dir <- Sys.getenv("FASTDER_FIG_DIR", file.path("results", "figures"))
out <- if (length(args) >= 1) args[[1]] else file.path(fig_dir, "fig_sim_jaccard.pdf")

source(file.path(dirname(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE))), "helpers.R"))
BASE_SIZE <- 8

jac <- read_result("config_full_simulation", "fuzzy_jaccard.csv")
fj <- pick_tool_param(jac, reference_pids(jac)) %>%
  mutate(jaccard = as.numeric(jaccard), scenario = relabel_scenario(scenario),
         event = relabel_event(sample),
         tool = factor(tool, levels = names(tool_labels)))

p <- ggplot(fj, aes(event, jaccard, fill = tool)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.25, width = 0.75,
               position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = tool_palette, labels = tool_labels, name = NULL) +
  facet_wrap(~ scenario, ncol = 1) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Event class", y = "Best exonic Jaccard per simulated transcript") +
  theme_pub()

ggsave(out, p, width = 7, height = 4.6)
cat("wrote", out, "\n")
