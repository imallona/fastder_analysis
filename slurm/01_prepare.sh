#!/bin/bash
# Cluster side of a split run: the timed configs up to their timed rules,
# then the accuracy-only configs in full. The timed rules run elsewhere.
#
#   sbatch slurm/01_prepare.sh
#
#SBATCH --job-name=fastder-prepare
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000
#SBATCH --output=slurm/logs/%x-%j.out
#SBATCH --signal=B:TERM@300

set -euo pipefail

repo="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
source "$repo/slurm/common.sh"

announce euler
echo
echo "=================== the plan ==================="
make euler "${MAKE_ARGS[@]}" EXTRA=-n 2>&1 | tail -30

echo
echo "=================== the run ==================="
run_targets euler

echo
echo "=================== rules over their declared memory ==================="
for benchmarks in workflow/logs/benchmarks/*/; do
    echo "$benchmarks"
    python3 workflow/scripts/check_declared_memory.py \
        --rules workflow/rules --benchmarks "$benchmarks" || true
done
