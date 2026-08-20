#!/bin/bash
# The simulation set at all four depths: 5M, 10M, 30M, 40M.
#
#   sbatch slurm/01_simulation.sh
#
# What this run covers, all of it new in the revision: the eight ASimulatoR
# event classes rather than four, chr21 and chr19 rather than chr21 alone, ten
# samples per scenario rather than five, the --no-stitch ablation, and every
# tool pinned to one core so the runtime comparison is like for like.
#
# The environments are built on a login node first, once. Compute nodes share a
# rate limited proxy, and 32 jobs solving environments at the same time abuse
# it:
#
#   module load eth_proxy
#   source /cluster/project/platt/$USER/miniforge3/bin/activate
#   make envs CONFIG=config_full_simulation.yaml EULER=1 \
#        CONDA_PREFIX_DIR=/cluster/project/platt/$USER/fastder-eval-envs
#
# This is the longest of the four run groups. run_asimulator asks for 24 hours
# per scenario and STAR indexes two chromosomes, so the driver's own wall clock
# is set well above the sum of what it waits for.
#
#SBATCH --job-name=fastder-sim
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000
#SBATCH --output=slurm/logs/%x-%j.out
#SBATCH --signal=B:TERM@300

source "$(dirname "$0")/common.sh"

announce simulations mjr-sweep
echo
echo "=================== the plan ==================="
make simulations mjr-sweep "${MAKE_ARGS[@]}" EXTRA=-n 2>&1 | tail -30

echo
echo "=================== the run ==================="
run_targets simulations
# The junction read-support sensitivity run reuses the simulated data of the
# 10M point, so it costs only its own fastder calls and has to follow them.
run_targets mjr-sweep
