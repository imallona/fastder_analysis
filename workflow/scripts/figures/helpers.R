# Theme, palette and panel builders shared by the assembled manuscript figures.
# Each assembly script sources this file and combines the panels with patchwork.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(patchwork)
})

BASE_SIZE <- 12

TOOLS <- c("fastder", "derfinder", "megadepth_baseline", "grohmm")

# Identical to the values in the fastder-evaluation reports, so figures match.
tool_palette <- c(
  derfinder          = "#66C2A5",
  fastder            = "#FC8D62",
  grohmm             = "#8DA0CB",
  megadepth_baseline = "#E78AC3"
)
tool_shapes <- c(fastder = 16, derfinder = 17, megadepth_baseline = 15, grohmm = 18)
# Distinct line types so near-coincident trend lines stay tellable apart.
tool_linetypes <- c(fastder = "solid", derfinder = "dashed",
                    megadepth_baseline = "dotted", grohmm = "dotdash")
# Small multiplicative offsets fan the tools out horizontally at each depth on
# the log axis, so near-identical values do not draw on top of each other.
tool_dodge <- c(derfinder = 0.92, fastder = 0.973,
                megadepth_baseline = 1.028, grohmm = 1.085)

# Strand classes use a grey/red/brown family, well outside the tools' pastel
# teal/orange/periwinkle/pink, so the two legends can never be confused.
strand_palette <- c(concordant = "#333333", discordant = "#d73027",
                    unstranded = "#bdbdbd", unmatched = "#8c510a")

tool_labels <- c(
  fastder = "fastder", derfinder = "derfinder",
  megadepth_baseline = "baseline", grohmm = "groHMM"
)

theme_pub <- function(base_size = BASE_SIZE) {
  theme_classic(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size),
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.subtitle = element_text(size = base_size - 1),
      plot.tag = element_text(face = "bold", size = base_size + 4)
    )
}

theme_pub_square <- function(base_size = BASE_SIZE) {
  theme_pub(base_size) + theme(aspect.ratio = 1)
}

# Dedicated legend grob: tools as colour + line type + point symbol (so it
# always matches the plotted lines), plus the strand-class fill. Built from a
# throwaway plot and extracted, because patchwork's guide collection rebuilds
# the keys as point-only and drops the line type.
make_combined_legend <- function() {
  suppressPackageStartupMessages(library(cowplot))
  tdf <- data.frame(tool = factor(names(tool_palette), levels = names(tool_palette)), x = 1, y = 1)
  cdf <- data.frame(category = factor(names(strand_palette), levels = names(strand_palette)), x = 1, y = 1)
  p <- ggplot() +
    geom_line(data = tdf, aes(x, y, colour = tool, linetype = tool)) +
    geom_point(data = tdf, aes(x, y, colour = tool, shape = tool), size = 2.6) +
    geom_tile(data = cdf, aes(x, y, fill = category)) +
    scale_colour_manual(values = tool_palette, labels = tool_labels, name = NULL) +
    scale_shape_manual(values = tool_shapes, labels = tool_labels, name = NULL) +
    scale_linetype_manual(values = tool_linetypes, labels = tool_labels, name = NULL) +
    scale_fill_manual(values = strand_palette, name = NULL) +
    theme_pub() + theme(legend.position = "bottom", legend.box = "vertical")
  # ggplot2 >= 3.5: get_legend returns an empty box; take the bottom guide-box.
  legend_grob <- tryCatch(
    cowplot::get_plot_component(p, "guide-box-bottom", return_all = TRUE),
    error = function(e) NULL)
  if (is.list(legend_grob) && !grid::is.grob(legend_grob)) {
    nonempty <- Filter(function(g) !inherits(g, "zeroGrob"), legend_grob)
    legend_grob <- if (length(nonempty) > 0) nonempty[[1]] else NULL
  }
  if (is.null(legend_grob) || inherits(legend_grob, "zeroGrob"))
    legend_grob <- cowplot::get_legend(p)
  legend_grob
}

# Parse mirrors benchmarks.Rmd: tool is encoded in the run_<tool> rule name.
load_benchmarks <- function(bench_dir) {
  files <- list.files(bench_dir, pattern = "\\.tsv$", recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) stop("no .tsv benchmark files under ", bench_dir)
  parse_one <- function(path) {
    rel <- sub(paste0(bench_dir, .Platform$file.sep), "", path)
    parts <- strsplit(rel, .Platform$file.sep, fixed = TRUE)[[1]]
    rule <- parts[1]
    d <- read_tsv(path, show_col_types = FALSE)
    d$rule <- rule
    d
  }
  df <- bind_rows(lapply(files, parse_one))
  df$wall_s <- as.numeric(df$s)
  df$max_rss_mb <- as.numeric(df$max_rss)
  df$tool <- dplyr::recode(df$rule,
    run_fastder = "fastder", run_derfinder = "derfinder",
    run_megadepth_baseline = "megadepth_baseline", run_grohmm = "grohmm",
    .default = NA_character_)
  df
}

# Large outlined point marks the per-tool median across invocations.
panel_speed <- function(bench_dir) {
  df <- load_benchmarks(bench_dir) %>% filter(tool %in% TOOLS)
  med <- df %>% group_by(tool) %>%
    summarise(wall_s = median(wall_s, na.rm = TRUE),
              max_rss_mb = median(max_rss_mb, na.rm = TRUE), .groups = "drop")
  ggplot(df, aes(wall_s, max_rss_mb, colour = tool)) +
    geom_point(size = 1.8, alpha = 0.35) +
    geom_point(data = med, aes(fill = tool), shape = 21, size = 3.4,
               colour = "black", stroke = 0.6, show.legend = FALSE) +
    scale_colour_manual(values = tool_palette, labels = tool_labels) +
    scale_fill_manual(values = tool_palette) +
    scale_x_log10() +
    scale_y_log10() +
    # No tool legend here; the shared legend comes from the line panels, whose
    # keys show the line type and point symbol together.
    guides(colour = "none") +
    labs(x = "Wall time per sample (s)", y = "Peak resident memory (MiB)") +
    theme_pub_square()
}

# Results tree to read. Defaults to the local mirror; the snakemake rule sets
# FASTDER_RESULTS_ROOT to the workflow results directory for a barbara rerun.
RESULTS_ROOT <- Sys.getenv("FASTDER_RESULTS_ROOT",
                           "/home/imallona/src/writing_fastder/barbara_results/results")
FIG_DIR <- Sys.getenv("FASTDER_FIG_DIR", "/home/imallona/src/writing_fastder/figures")

# Read one CSV from every config_full_simulation* run, stamping the depth its
# name encodes (the base config is 10M). Mirrors the loader in meta.Rmd.
load_depth_sweep <- function(file_name, root = RESULTS_ROOT) {
  run_dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  run_dirs <- run_dirs[grepl("^config_full_simulation", basename(run_dirs))]
  depth_of <- function(d) {
    m <- str_match(d, "_([0-9]+)M$")[, 2]
    if (is.na(m)) 10L else as.integer(m)
  }
  read_run <- function(dir_path) {
    path <- file.path(dir_path, file_name)
    if (!file.exists(path)) return(NULL)
    read_csv(path, show_col_types = FALSE) %>% mutate(depth_M = depth_of(basename(dir_path)))
  }
  bind_rows(lapply(run_dirs, read_run))
}

# Panel: gffcompare sensitivity and precision against depth, exon and
# transcript levels, averaged over samples and parameters. Verbatim from the
# meta.Rmd sens_prec chunk, restyled to the shared clean theme.
panel_depth <- function(which_levels = c("Transcript", "Exon")) {
  summary_all <- load_depth_sweep("summary.csv")
  levels_long <- bind_rows(
    summary_all %>% transmute(tool, scenario, depth_M, level = "Transcript",
                              sensitivity = transcript_sens, precision = transcript_prec),
    summary_all %>% transmute(tool, scenario, depth_M, level = "Exon",
                              sensitivity = exon_sens, precision = exon_prec)
  ) %>%
    filter(level %in% which_levels) %>%
    pivot_longer(c(sensitivity, precision), names_to = "metric", values_to = "value") %>%
    group_by(tool, scenario, depth_M, level, metric) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(scenario = relabel_scenario(scenario),
           depth_x = depth_M * tool_dodge[as.character(tool)])
  ggplot(levels_long, aes(depth_x, value, colour = tool, shape = tool, linetype = tool)) +
    geom_line(linewidth = 0.7) + geom_point(size = 2.4) +
    scale_x_continuous(trans = "log10", breaks = sort(unique(levels_long$depth_M))) +
    scale_colour_manual(values = tool_palette, labels = tool_labels) +
    scale_shape_manual(values = tool_shapes, labels = tool_labels) +
    scale_linetype_manual(values = tool_linetypes, labels = tool_labels) +
    coord_cartesian(ylim = c(0, 100)) +
    # Single level (the main figure uses exon only): put the short metric names
    # on top and the long scenario label on the rotated right strip, so the
    # square facets do not clip. Both levels: keep the report's layout.
    (if (length(which_levels) == 1)
       facet_grid(scenario ~ metric, labeller = labeller(
         scenario = label_wrap_gen(12),
         metric = c(precision = "Precision", sensitivity = "Sensitivity")))
     else
       facet_grid(level + metric ~ scenario, labeller = labeller(
         scenario = label_wrap_gen(12),
         metric = c(precision = "Precision", sensitivity = "Sensitivity")))) +
    labs(x = "Reads per sample (M)", y = "Percent") +
    theme_pub_square() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(colour = "none", shape = "none", linetype = "none")
}

# Panel: fraction of exon boundaries within 5 bp of a reference boundary,
# against depth. Verbatim from the meta.Rmd boundary chunk, restyled.
panel_boundary <- function() {
  distances_all <- load_depth_sweep("fuzzy_distances.csv")
  b5 <- distances_all %>%
    mutate(distance = as.integer(distance)) %>%
    group_by(tool, scenario, depth_M, sample, param_id) %>%
    summarise(pct = mean(abs(distance) <= 5) * 100, .groups = "drop") %>%
    group_by(tool, scenario, depth_M) %>%
    summarise(pct = mean(pct), .groups = "drop") %>%
    mutate(scenario = relabel_scenario(scenario),
           depth_x = depth_M * tool_dodge[as.character(tool)])
  ggplot(b5, aes(depth_x, pct, colour = tool, shape = tool, linetype = tool)) +
    geom_line(linewidth = 0.7) + geom_point(size = 2.4) +
    scale_x_continuous(trans = "log10", breaks = sort(unique(b5$depth_M))) +
    scale_colour_manual(values = tool_palette, labels = tool_labels) +
    scale_shape_manual(values = tool_shapes, labels = tool_labels) +
    scale_linetype_manual(values = tool_linetypes, labels = tool_labels) +
    coord_cartesian(ylim = c(0, 100)) +
    facet_wrap(~ scenario, labeller = label_wrap_gen(12)) +
    labs(x = "Reads per sample (M)", y = "Exon boundaries within\n5 bp of truth (%)") +
    theme_pub_square() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(colour = "none", shape = "none", linetype = "none")
}

# --- summary.Rmd helpers, ported verbatim so strand and multi-exon panels
# select the same parameter combinations the report figures used. ---

scenario_labels <- c(template_and_variant = "Reference and variant", variant_only = "Variant")
relabel_scenario <- function(v) {
  factor(v, levels = names(scenario_labels), labels = scenario_labels[names(scenario_labels)])
}

extract_mc <- function(pid) suppressWarnings(as.numeric(str_extract(pid, "(?<=mc)[0-9.]+")))
extract_pt <- function(pid) {
  v <- suppressWarnings(as.numeric(str_extract(pid, "(?<=pt)[0-9.]+")))
  ifelse(is.na(v), 0, v)
}

# Best parameter per tool by median Jaccard, matching baselines to fastder on
# the shared mc (and pt for derfinder) axes; grohmm at its own best.
best_pids <- function(jaccard) {
  best_fastder <- jaccard %>% filter(tool == "fastder") %>%
    mutate(jaccard = as.numeric(jaccard)) %>%
    group_by(param_id) %>% summarise(m = median(jaccard), .groups = "drop") %>%
    slice_max(m, n = 1, with_ties = FALSE) %>% pull(param_id)
  mc <- extract_mc(best_fastder); pt <- extract_pt(best_fastder)
  best_grohmm <- jaccard %>% filter(tool == "grohmm") %>%
    mutate(jaccard = as.numeric(jaccard)) %>%
    group_by(param_id) %>% summarise(m = median(jaccard), .groups = "drop") %>%
    slice_max(m, n = 1, with_ties = FALSE) %>% pull(param_id)
  if (length(best_grohmm) == 0) best_grohmm <- NA_character_
  list(fastder = best_fastder,
       megadepth_baseline = sprintf("mc%s", format(mc, trim = TRUE)),
       derfinder = sprintf("mc%s_pt%s", format(mc, trim = TRUE), format(pt, trim = TRUE)),
       grohmm = best_grohmm)
}

pick_tool_param <- function(d, pids) {
  d %>% filter((tool == "fastder" & param_id == pids$fastder) |
               (tool == "megadepth_baseline" & param_id == pids$megadepth_baseline) |
               (tool == "derfinder" & param_id == pids$derfinder) |
               (tool == "grohmm" & param_id == pids$grohmm))
}

read_result <- function(config, file_name) {
  read_csv(file.path(RESULTS_ROOT, config, file_name), show_col_types = FALSE)
}

# Panel: per-tool strand concordance. Verbatim from summary.Rmd tool_strand,
# restyled. Only fastder assigns a strand; the others are all unstranded.
# compact pools the five AS-event samples into two scenario cells (main figure);
# the default keeps the full per-sample grid (supplement).
panel_strand <- function(config = "config_full_simulation", compact = FALSE) {
  pids <- best_pids(read_result(config, "fuzzy_jaccard.csv"))
  fs <- read_result(config, "fuzzy_strand.csv") %>%
    pick_tool_param(pids) %>%
    mutate(n = as.integer(n_fastder_transcripts), scenario = relabel_scenario(scenario))
  if (compact) {
    fs_pct <- fs %>% group_by(scenario, tool, category) %>%
      summarise(n = sum(n), .groups = "drop") %>%
      group_by(scenario, tool) %>% mutate(pct = n / sum(n) * 100) %>% ungroup()
  } else {
    fs_pct <- fs %>% group_by(scenario, sample, tool) %>%
      mutate(pct = n / sum(n) * 100) %>% ungroup()
  }
  fs_pct$category <- factor(fs_pct$category,
                            levels = c("concordant", "discordant", "unstranded", "unmatched"))
  fs_pct$tool <- factor(fs_pct$tool, levels = rev(TOOLS))
  p <- ggplot(fs_pct, aes(pct, tool, fill = category)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = strand_palette) +
    scale_y_discrete(labels = tool_labels) +
    labs(x = "Fraction of tool transcripts (%)", y = NULL) +
    theme_pub() + theme(legend.position = "none")
  if (compact)
    p + facet_wrap(~ scenario, nrow = 1, labeller = label_wrap_gen(20))
  else
    p + facet_grid(scenario ~ sample, labeller = labeller(scenario = label_wrap_gen(12)))
}

# Panel: multi-exon transcripts produced per tool. Verbatim from summary.Rmd
# tool_multiexon, restyled. Baselines emit none, pinned to 1 on the log axis.
# No per-point labels (they collided); tools read from the shared legend.
panel_multiexon <- function(config = "config_full_simulation") {
  pids <- best_pids(read_result(config, "fuzzy_jaccard.csv"))
  d <- read_result(config, "summary.csv") %>%
    pick_tool_param(pids) %>%
    distinct(tool, scenario, param_id, query_mrnas, query_multi_exon) %>%
    mutate(query_mrnas = as.integer(query_mrnas),
           query_multi_exon = as.integer(query_multi_exon),
           scenario = relabel_scenario(scenario)) %>%
    group_by(tool, scenario) %>%
    summarise(query_mrnas = sum(query_mrnas),
              query_multi_exon = sum(query_multi_exon), .groups = "drop") %>%
    mutate(query_multi_exon_plot = pmax(query_multi_exon, 1L))
  ggplot(d, aes(query_mrnas, query_multi_exon_plot, colour = tool, shape = tool)) +
    geom_point(size = 3.2, alpha = 0.9) +
    scale_colour_manual(values = tool_palette, labels = tool_labels) +
    scale_shape_manual(values = tool_shapes, labels = tool_labels) +
    scale_x_log10() + scale_y_log10() +
    facet_wrap(~ scenario, ncol = 2, labeller = label_wrap_gen(12)) +
    labs(x = "Predicted transcripts (count)", y = "Multi-exon transcripts (count)") +
    theme_pub_square() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
}

# Panel: cumulative fraction of predicted boundaries within each distance of a
# reference boundary. Verbatim from summary.Rmd tool_boundary_ecdf, restyled;
# a steep left rise means boundaries land on splice sites.
panel_boundary_cdf <- function(config = "config_full_simulation") {
  pids <- best_pids(read_result(config, "fuzzy_jaccard.csv"))
  fd <- read_result(config, "fuzzy_distances.csv") %>%
    pick_tool_param(pids) %>%
    mutate(abs_distance = abs(as.integer(distance)),
           scenario = relabel_scenario(scenario),
           tool = factor(tool, levels = c("megadepth_baseline", "derfinder", "grohmm", "fastder")))
  ggplot(fd, aes(abs_distance, colour = tool, linetype = tool)) +
    stat_ecdf(geom = "step", linewidth = 0.8, alpha = 0.85) +
    scale_colour_manual(values = tool_palette, labels = tool_labels) +
    scale_linetype_manual(values = tool_linetypes, labels = tool_labels) +
    scale_x_continuous(trans = "log1p", breaks = c(0, 5, 100, 10000),
                       labels = c("0", "5", "100", "10k")) +
    coord_cartesian(ylim = c(0, 1)) +
    facet_wrap(~ scenario, labeller = label_wrap_gen(12)) +
    labs(x = "Distance to nearest true boundary (bp)", y = "Cumulative fraction of boundaries") +
    theme_pub_square() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(colour = "none", shape = "none", linetype = "none")
}

# --- Figure 2: GTEx chr19 tool-comparison panels, re-plotted from the local
# config_gtex_comparison CSVs (summary.csv and the fastder chain_stats.csv). ---

GTEXCMP <- "config_gtex_comparison"
FASTDER_ROSE <- unname(tool_palette["fastder"])

# Tissue colours, chosen clear of the tool palette (teal/orange/periwinkle/pink).
tissue_palette <- c(brain = "#6a51a3", heart = "#cb181d",
                    muscle = "#41ab5d", blood = "#2171b5")

# Panel: novel expressed-region exons per tissue (no overlapping Ensembl exon),
# genome-wide GTEx. Counts are the gtex_concordance run output. Portrait bars.
NOVEL_CSV <- file.path(FIG_DIR, "novel_exons.csv")
panel_novel <- function(csv = NOVEL_CSV) {
  d <- read_csv(csv, show_col_types = FALSE) %>%
    mutate(tissue = factor(tissue, levels = c("brain", "heart", "muscle", "blood")),
           lab = paste0(novel, "\n(", pct, "%)"))
  ggplot(d, aes(tissue, novel, fill = tissue)) +
    geom_col(width = 0.8) +
    geom_text(aes(label = lab), vjust = -0.2, size = 3, lineheight = 0.85) +
    scale_fill_manual(values = tissue_palette, guide = "none") +
    expand_limits(y = max(d$novel) * 1.25) +
    labs(x = NULL, y = "Novel ER exons (count)") +
    theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# Median precision per tool across the 32 sub-groups (chr19). `level` picks the
# gffcompare level: "exon" or "transcript".
panel_gtexcmp_precision <- function(config = GTEXCMP, level = "exon") {
  col <- if (level == "transcript") "transcript_prec" else "exon_prec"
  d <- read_result(config, "summary.csv") %>%
    group_by(tool) %>% summarise(prec = median(.data[[col]], na.rm = TRUE), .groups = "drop") %>%
    filter(tool %in% TOOLS)
  ggplot(d, aes(reorder(tool, -prec), prec, fill = tool)) +
    geom_col(width = 0.75) +
    scale_fill_manual(values = tool_palette, guide = "none") +
    scale_x_discrete(labels = tool_labels) +
    labs(x = NULL, y = paste0(tools::toTitleCase(level), " precision (%)")) +
    theme_pub() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# Exons per stitched fastder region (chr19), capped at 6+.
panel_gtexcmp_exons <- function(config = GTEXCMP) {
  d <- read_result(config, "chain_stats.csv") %>%
    mutate(n_exons = pmin(as.integer(n_exons), 6)) %>% count(n_exons) %>%
    mutate(lab = ifelse(n_exons == 6, "6+", as.character(n_exons)))
  ggplot(d, aes(factor(n_exons), n)) +
    geom_col(fill = FASTDER_ROSE, width = 0.85) +
    scale_x_discrete(labels = d$lab) + scale_y_log10() +
    labs(x = "Exons per region\n(count)", y = "Number of regions") +
    theme_pub()
}

# Exonic length per fastder region (chr19).
panel_gtexcmp_length <- function(config = GTEXCMP) {
  d <- read_result(config, "chain_stats.csv") %>%
    mutate(total_exon_length = as.numeric(total_exon_length))
  ggplot(d, aes(total_exon_length)) +
    geom_histogram(bins = 40, fill = FASTDER_ROSE) +
    scale_x_log10() +
    labs(x = "Exonic length per region (bp)", y = "Number of regions") +
    theme_pub()
}

# Stitched-expressed-region score per fastder region (chr19).
panel_gtexcmp_score <- function(config = GTEXCMP) {
  d <- read_result(config, "chain_stats.csv") %>% mutate(score = as.numeric(score))
  ggplot(d, aes(score)) +
    geom_histogram(bins = 40, fill = FASTDER_ROSE) +
    scale_x_log10() +
    labs(x = "Expressed-region score", y = "Number of regions") +
    theme_pub()
}

# fastder regions per chromosome and strand, genome-wide (concordance run spans
# all chromosomes, the chr19 comparison run does not). Sub-groups pooled per tissue.
GTEXCONC <- "config_gtex_concordance"
genomic_strand_fill <- c("-" = "#F8766D", "." = "#00BA38", "+" = "#619CFF")
panel_gtex_genomic_dist <- function(config = GTEXCONC) {
  chrom_order <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")
  d <- read_result(config, "chain_stats.csv") %>%
    mutate(tissue = factor(sub("_[0-9]+$", "", scenario),
                           levels = c("brain", "heart", "muscle", "blood")),
           chrom = factor(chrom, levels = chrom_order),
           strand = factor(strand, levels = c("-", ".", "+"))) %>%
    filter(!is.na(chrom))
  ggplot(d, aes(chrom, fill = strand)) +
    geom_bar() +
    scale_fill_manual(values = genomic_strand_fill, name = "strand") +
    facet_wrap(~ tissue, ncol = 1, scales = "free_y") +
    labs(x = NULL, y = "fastder regions (count)") +
    theme_pub() +
    theme(legend.position = "right", legend.title = element_text(),
          axis.text.x = element_text(angle = 45, hjust = 1))
}

# Panel: troponin marker loci as compact per-tissue ER tracks. One facet per
# gene, one row per tissue, ER exons drawn as boxes. Tissue-restricted: TNNT2
# and TNNI3 in heart, TNNT3 in skeletal muscle. Data extracted from the
# per-sub-group GTFs into marker_loci.csv.
MARKER_CSV <- file.path(FIG_DIR, "marker_loci.csv")
panel_marker_loci <- function(csv = MARKER_CSV) {
  tissues <- c("brain", "heart", "muscle", "blood")
  d <- read_csv(csv, show_col_types = FALSE) %>%
    mutate(tissue = factor(tissue, levels = rev(tissues)),
           gene = factor(gene, levels = c("TNNT2", "TNNT3", "TNNI3")),
           start = start / 1e6, end = end / 1e6)
  ggplot(d, aes(xmin = start, xmax = end,
                ymin = as.integer(tissue) - 0.4, ymax = as.integer(tissue) + 0.4)) +
    geom_rect(fill = FASTDER_ROSE) +
    scale_y_continuous(breaks = seq_along(tissues), labels = rev(tissues),
                       limits = c(0.5, length(tissues) + 0.5)) +
    facet_wrap(~ gene, scales = "free_x", nrow = 1, labeller = labeller(
      gene = c(TNNT2 = "TNNT2 (chr1)", TNNT3 = "TNNT3 (chr11)", TNNI3 = "TNNI3 (chr19)"))) +
    labs(x = "position (Mb)", y = NULL) +
    theme_pub()
}

# Embed a schematic PDF as a panel, rasterised at 600 dpi. rasterGrob letterboxes
# to the cell; a grImport2 vector embed was tried but overflows narrow cells.
# width_fill stretches to full panel width instead of centring at native aspect.
wrap_pdf <- function(pdf_path, dpi = 600, width_fill = FALSE) {
  tmp <- tempfile(fileext = "")
  system2("pdftoppm", c("-png", "-r", dpi, "-singlefile", pdf_path, tmp))
  img <- png::readPNG(paste0(tmp, ".png"))
  grob <- if (width_fill) {
    grid::rasterGrob(img, width = grid::unit(1, "npc"), interpolate = TRUE)
  } else {
    grid::rasterGrob(img, interpolate = TRUE)
  }
  wrap_elements(full = grob)
}

wrap_png <- function(name, fig_dir = FIG_DIR) {
  wrap_elements(full = grid::rasterGrob(
    png::readPNG(file.path(fig_dir, name)), interpolate = TRUE))
}

panel_placeholder <- function(text) {
  ggplot() + annotate("text", 0, 0, label = text, size = 4) +
    theme_void()
}
