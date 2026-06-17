# Per-event Jaccard over the parameter grid. min_length on x; min_coverage and
# position_tolerance as facet rows so pt invariance shows as near-identical bands.
# Reads config_full_simulation; rebuilds from the local mirror.

args <- commandArgs(trailingOnly = TRUE)
fig_dir <- Sys.getenv("FASTDER_FIG_DIR", "/home/imallona/src/writing_fastder/figures")
out <- if (length(args) >= 1) args[[1]] else file.path(fig_dir, "fig_sim_event_jaccard.png")

source(file.path(dirname(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE))), "helpers.R"))

as_event_levels <- c("es", "mes", "afe", "ale", "mixed")
as_event_labels <- c(es = "exon skipping (ES)", mes = "multiple exon skipping (MES)",
  afe = "alternative first exon (AFE)", ale = "alternative last exon (ALE)", mixed = "mixed events")
label_as <- function(v) factor(v, levels = as_event_levels, labels = as_event_labels[as_event_levels])
scenario_palette <- c("Reference and variant" = "#1b9e77", "Variant" = "#d95f02")
parse_param_id <- function(d) d %>% mutate(
  min_coverage = as.numeric(str_extract(param_id, "(?<=mc)[0-9.]+")),
  min_length = as.numeric(str_extract(param_id, "(?<=ml)[0-9.]+")),
  position_tolerance = { v <- as.numeric(str_extract(param_id, "(?<=pt)[0-9.]+")); ifelse(is.na(v), 0, v) })

fj <- read_result("config_full_simulation", "fuzzy_jaccard.csv") %>%
  filter(tool == "fastder") %>%
  mutate(jaccard = as.numeric(jaccard), as_event = label_as(sample),
         scenario = relabel_scenario(scenario)) %>%
  parse_param_id() %>% mutate(param_sub = paste0("ml", min_length))

p <- ggplot(fj, aes(param_sub, jaccard, fill = scenario)) +
  geom_boxplot(outlier.size = 0.25, linewidth = 0.3) +
  scale_fill_manual(values = scenario_palette, name = NULL) +
  labs(x = "min_length", y = "Best exonic Jaccard against reference") +
  coord_cartesian(ylim = c(0, 1)) +
  facet_grid(min_coverage + position_tolerance ~ as_event,
             labeller = labeller(min_coverage = function(v) paste0("mc ", v),
                                 position_tolerance = function(v) paste0("pt ", v))) +
  theme_pub() + theme(legend.position = "bottom")

ggsave(out, p, width = 12, height = 11, dpi = 200)
cat("wrote", out, "\n")
