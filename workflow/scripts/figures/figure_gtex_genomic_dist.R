# Genome-wide genomic distribution of fastder regions per chromosome and strand.
# Reads config_gtex_concordance/chain_stats.csv.

args <- commandArgs(trailingOnly = TRUE)
fig_dir <- Sys.getenv("FASTDER_FIG_DIR", "/home/imallona/src/writing_fastder/figures")
out <- if (length(args) >= 1) args[[1]] else file.path(fig_dir, "fig_gtexcmp_genomic_dist.pdf")

source(file.path(dirname(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE))), "helpers.R"))

p <- panel_gtex_genomic_dist()
ggsave(out, p, width = 9, height = 9)
ggsave(sub("\\.pdf$", ".png", out), p, width = 9, height = 9, dpi = 300)
cat("wrote", out, "and", sub("\\.pdf$", ".png", out), "\n")
