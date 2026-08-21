#!/bin/bash
# The two TDP-43 runs, showcase and panel.
#
#   sbatch slurm/03_tdp43.sh
#
# Cheap: four recount3 BigWigs over chr8, chr19 and chr20. Rerun because
# whole-file library size changes what a CPM threshold means.
#
# The configs still carry the old hand-tuned thresholds, 1.0 and 0.02 CPM.
# This run feeds the threshold sweep; it is not a finished result.
#
#SBATCH --job-name=fastder-tdp43
#SBATCH --time=24:00:00
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

announce tdp43 tdp43-panel
echo
echo "=================== the plan ==================="
make tdp43 tdp43-panel "${MAKE_ARGS[@]}" EXTRA=-n 2>&1 | tail -30

echo
echo "=================== the run ==================="
run_targets tdp43
run_targets tdp43-panel
