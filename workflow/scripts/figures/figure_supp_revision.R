# Supplementary figures added in the revision, one file each so the manuscript
# can place them independently:
#
#   supp_ablation.pdf              junction integration on and off, against depth
#   supp_min_junction_reads.pdf    accuracy against the junction read-support filter
#   supp_scaling.pdf               wall time and peak memory against cores
#
# Each panel reads the tidy CSV its collector script wrote and saves the frame
# it drew next to the figure, so every plotted number is on disk.
#
# Usage: Rscript figure_supp_revision.R <out_dir>

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else "."

source(file.path(dirname(sub("--file=", "",
  grep("--file=", commandArgs(FALSE), value = TRUE))), "helpers.R"))

# Structural ceilings of the scaling workload, annotated on the panel. The
# simulation loads ten samples over two chromosomes; override for another
# workload with the two environment variables.
scaling_samples <- as.integer(Sys.getenv("FASTDER_SCALING_SAMPLES", "10"))
scaling_chroms <- as.integer(Sys.getenv("FASTDER_SCALING_CHROMS", "2"))

save_both <- function(plot, name, width, height) {
  pdf_path <- file.path(out_dir, paste0(name, ".pdf"))
  ggsave(pdf_path, plot, width = width, height = height, limitsize = FALSE)
  ggsave(sub("\\.pdf$", ".svg", pdf_path), plot, width = width, height = height,
         limitsize = FALSE)
  cat("wrote", pdf_path, "\n")
}

save_both(panel_ablation(file.path(out_dir, "ablation.csv")),
          "supp_ablation", width = 7.0, height = 6.0)

save_both(panel_min_junction_reads(file.path(out_dir, "min_junction_reads.csv")),
          "supp_min_junction_reads", width = 7.0, height = 3.6)

save_both(panel_scaling(file.path(out_dir, "scaling.csv"),
                        samples = scaling_samples, chromosomes = scaling_chroms),
          "supp_scaling", width = 7.0, height = 3.6)
