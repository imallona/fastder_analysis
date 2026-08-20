#!/bin/bash
# Step 0 on Euler: check the cluster behaves the way the run scripts assume.
# Costs one short job. Run it once per cluster, and again after a Slurm or
# executor plugin upgrade.
#
#   sbatch slurm/00_probe.sh
#
# It answers four things a login node cannot:
#
# 1. Does the environment come up inside a batch job: conda activates, the
#    Snakefile parses, the DAG builds.
# 2. May a batch job submit to Slurm from here. If not, the driver has to sit
#    on a login node and everything below changes.
# 3. Does the benchmark directive still produce a TSV when the payload runs on
#    a compute node. Every runtime number in the paper depends on it.
# 4. Are the EPYC_7763 nodes the profile pins the timed rules to reachable
#    from this account.
#
#SBATCH --job-name=fastder-probe
#SBATCH --time=00:45:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2000
#SBATCH --output=slurm/logs/%x-%j.out

set -euo pipefail

module load eth_proxy

CONDA_INIT="${CONDA_INIT:-/cluster/project/platt/$USER/miniforge3/bin/activate}"
CONDA_ENV="${CONDA_ENV:-snakemake}"
CONDA_PREFIX_DIR="${CONDA_PREFIX_DIR:-/cluster/project/platt/$USER/fastder-eval-envs}"

source "$CONDA_INIT"
conda activate "$CONDA_ENV"

cd "${SLURM_SUBMIT_DIR:-$PWD}"

echo "=================== the driver's own environment ==================="
echo "host      $(hostname)"
echo "job id    ${SLURM_JOB_ID:-none}"
echo "snakemake $(snakemake --version)"
python -c 'import snakemake_executor_plugin_slurm as m; print("slurm plugin", m.__version__)'
echo

echo "=================== 1. can it plan from here ==================="
make smoke EULER=1 CONDA_INIT="$CONDA_INIT" CONDA_ENV="$CONDA_ENV" \
     CONDA_PREFIX_DIR="$CONDA_PREFIX_DIR" EXTRA=-n 2>&1 | tail -25 || true
echo

echo "=================== 2. and 3. can it submit, and does benchmarking survive ==================="
# One real child job that writes a benchmark TSV. record_host_info is the
# cheapest rule that produces one, and its output is what Methods reports as
# the benchmark machine.
make smoke EULER=1 CONDA_INIT="$CONDA_INIT" CONDA_ENV="$CONDA_ENV" \
     CONDA_PREFIX_DIR="$CONDA_PREFIX_DIR" \
     EXTRA="--until record_host_info"
echo
echo "host_info.tsv:"
cat workflow/results/config_quick_light/host_info.tsv || echo "not written"
echo

echo "=================== 4. are the pinned nodes reachable ==================="
sinfo -h -o "%D %P %t" --constraint=EPYC_7763 -p normal.4h,normal.24h | head
echo
echo "If nothing is listed above, the --constraint in profiles/euler/config.yaml"
echo "names a model this account cannot reach and the timed rules will queue"
echo "forever. Pick another model from sinfo -o '%f' and change it in one place."
