#!/bin/bash
# The simulation set at all four depths: 5M, 10M, 30M, 40M.
#
#   sbatch slurm/01_simulation.sh
#
# New here: eight event classes, chr21 and chr19, ten samples per scenario,
# the --no-stitch ablation, every tool on one core.
#
# Build the environments once on a login node first; the proxy is shared:
#
#   module load eth_proxy
#   source /cluster/project/platt/$USER/miniforge3/bin/activate
#   make envs CONFIG=config_full_simulation.yaml EULER=1 \
#        CONDA_PREFIX_DIR=/cluster/project/platt/$USER/fastder-eval-envs
#
# The longest of the four run groups.
#
#SBATCH --job-name=fastder-sim
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

announce simulations mjr-sweep
echo
echo "=================== the plan ==================="
make simulations mjr-sweep "${MAKE_ARGS[@]}" EXTRA=-n 2>&1 | tail -30

echo
echo "=================== the run ==================="
run_targets simulations
# Reuses the 10M simulated data, so it has to follow.
run_targets mjr-sweep
