#!/bin/bash
# The two GTEx runs: the chr19 four-tool comparison, then the genome-wide
# atlas with the core scaling sweep.
#
#   sbatch slurm/02_gtex.sh
#
# Both are rerun because library size is now whole-file, which changes what a
# CPM threshold means. chr19 gives the single-core cross-tool runtimes.
#
# Downloads dominate the first hours: 160 BigWigs at 124 MB, plus one junction
# matrix per tissue at 2.8 GB gzipped, through the shared proxy.
#
# About 100 GB, needed for weeks, so project storage rather than scratch.
#
#SBATCH --job-name=fastder-gtex
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000
#SBATCH --output=slurm/logs/%x-%j.out
#SBATCH --signal=B:TERM@300

set -euo pipefail

# Slurm copies the batch script into a spool directory before running it, so
# $0 does not point into the repo. sbatch is called from the repo root, which
# is what SLURM_SUBMIT_DIR holds.
repo="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
source "$repo/slurm/common.sh"

announce gtex-comparison gtex
echo
echo "=================== the plan ==================="
make gtex-comparison gtex "${MAKE_ARGS[@]}" EXTRA=-n 2>&1 | tail -30

echo
echo "=================== the run ==================="
# The comparison first: the runtime table depends on it.
run_targets gtex-comparison
run_targets gtex
