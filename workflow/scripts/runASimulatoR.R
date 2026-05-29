# Redirect all output to Snakemake log file
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# ASimulatoR reads the GTF and FASTAs from a single directory and caches an
# exon-superset .rda next to the GTF. The reference inputs are passed as
# explicit files; symlink them into a private temp directory and point
# ASimulatoR there, so the .rda cache lands in scratch and never in the
# shared, scope-keyed reference folder. A write into the reference folder
# would change file timestamps there and re-trigger every downstream rule.
gtf <- snakemake@input[["gtf"]]
fastas <- snakemake@input[["fastas"]]
outdir <- snakemake@params[["outdir"]]

reference_dir <- dirname(normalizePath(gtf))
input_dir <- tempfile("asim_ref_")
dir.create(input_dir)
file.symlink(normalizePath(gtf), file.path(input_dir, basename(gtf)))
for (fa in fastas) {
  file.symlink(normalizePath(fa), file.path(input_dir, basename(fa)))
}

# Simulation parameters (from config.yaml via Snakefile params)
seed <- snakemake@params[["seed"]]
ncores <- snakemake@threads
seq_depth <- snakemake@params[["seq_depth"]]
multi_events_per_exon <- snakemake@params[["multi_events_per_exon"]]
strand_specific <- snakemake@params[["strand_specific"]]
probs_as_freq <- snakemake@params[["probs_as_freq"]]

# Python dict (config.yaml) -> R named list -> named numeric vector
event_probs <- unlist(snakemake@params[["events"]])

# Run simulation
library(ASimulatoR)

simulate_alternative_splicing(
  input_dir = input_dir,
  outdir = outdir,
  event_probs = event_probs,
  num_reps = c(1, 0),
  seed = seed,
  ncores = ncores,
  seq_depth = seq_depth,
  multi_events_per_exon = multi_events_per_exon,
  strand_specific = strand_specific,
  probs_as_freq = probs_as_freq
)

# The staged reference is scratch; drop it once the simulation has succeeded.
unlink(input_dir, recursive = TRUE)

# Save metadata (written after successful simulation only)
library(yaml)

metadata <- list(
  timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  params = list(
    reference_dir = reference_dir,
    seed = seed,
    ncores = ncores,
    seq_depth = seq_depth,
    multi_events_per_exon = multi_events_per_exon,
    strand_specific = strand_specific,
    probs_as_freq = probs_as_freq,
    event_probs = as.list(event_probs)
  )
)

write_yaml(metadata, file.path(outdir, "simulation_metadata.yaml"))
