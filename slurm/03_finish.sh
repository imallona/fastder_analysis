#!/bin/bash
# Evaluation and reports of the timed configs, one Slurm job per rule, then
# the figures. Runs after 02_timed.sh.
#
#   sbatch slurm/03_finish.sh
#
# The cross-depth report is rendered in this job.
#
#SBATCH --job-name=fastder-finish
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4000
#SBATCH --output=slurm/logs/%x-%j.out
#SBATCH --signal=B:TERM@300

set -euo pipefail

repo="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
source "$repo/slurm/common.sh"

announce timed-configs PASSES=rest
echo
echo "=================== the plan ==================="
make timed-configs PASSES=rest "${MAKE_ARGS[@]}" EXTRA=-n 2>&1 | tail -30

echo
echo "=================== the run ==================="
run_targets timed-configs PASSES=rest
run_targets figures
