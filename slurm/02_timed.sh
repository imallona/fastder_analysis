#!/bin/bash
# Timed rules of the prepared configs, one job at a time in this allocation.
#
#   sbatch slurm/02_timed.sh                      # every timed config
#   sbatch slurm/02_timed.sh tdp43 tdp43-panel    # some of them
#
# The node is shared with other jobs; host_info.tsv has its load.
# --cpus-per-task is the largest scaling_cores, --constraint the CPU model of
# the timed rules in the profile.
#
#SBATCH --job-name=fastder-timed
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4500
#SBATCH --constraint=EPYC_7763
#SBATCH --output=slurm/logs/%x-%j.out
#SBATCH --signal=B:TERM@300

set -euo pipefail

# Rules run in this allocation.
EULER=

repo="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
source "$repo/slurm/common.sh"

targets=("${@:-timed-configs}")
cores="${SLURM_CPUS_PER_TASK:-16}"
MAKE_ARGS+=(
    PASSES="check timed"
    QUIET_LOAD=
    CORES="$cores"
    TIMED_CORES="$cores"
    MEM_MB="$((cores * ${SLURM_MEM_PER_CPU:-4500}))"
)

announce "${targets[@]}"
echo
run_targets "${targets[@]}"
