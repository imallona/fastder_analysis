#!/bin/bash
# Evaluation and reports of the timed configs, one Slurm job per rule, then
# the figures. Runs after 02_timed.sh and stops at a config with a timed
# rule left to run.
#
#   sbatch slurm/03_finish.sh                      # every timed config, figures
#   sbatch slurm/03_finish.sh tdp43 tdp43-panel    # some of them, no figures
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

targets=("${@:-timed-configs}")
MAKE_ARGS+=(PASSES="untimed rest")

announce "${targets[@]}"
echo
run_targets "${targets[@]}"
if [ "$#" -eq 0 ]; then
    run_targets figures
fi
