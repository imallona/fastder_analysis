input_dir <- snakemake@input[["reference"]]
outdir    <- snakemake@params[["outdir"]]

# Simulation parameters (from config.yaml via Snakefile params)
seed                  <- snakemake@params[["seed"]]
ncores                <- snakemake@params[["ncores"]]
seq_depth             <- snakemake@params[["seq_depth"]]
multi_events_per_exon <- snakemake@params[["multi_events_per_exon"]]
strand_specific       <- snakemake@params[["strand_specific"]]
probs_as_freq         <- snakemake@params[["probs_as_freq"]]

# Python dict (config.yaml) -> R named list -> named numeric vector
event_probs <- unlist(snakemake@params[["events"]])

# Run simulation
library(ASimulatoR)

simulate_alternative_splicing(
  input_dir             = input_dir,
  outdir                = outdir,
  event_probs           = event_probs,
  num_reps              = c(1, 0),  # 1 sample per run, no second group
  seed                  = seed,
  ncores                = ncores,
  seq_depth             = seq_depth,
  multi_events_per_exon = multi_events_per_exon,
  strand_specific       = strand_specific,
  probs_as_freq         = probs_as_freq
)
