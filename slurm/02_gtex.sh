#!/bin/bash
# The two GTEx runs: the four-tool comparison on chr19, then the genome-wide
# atlas, which also carries the core scaling sweep.
#
#   sbatch slurm/02_gtex.sh
#
# Both are rerun because library size is now the whole file rather than the
# analysed chromosomes, which changes what a CPM threshold means. The chr19 run
# is also where the single-core cross-tool runtimes come from.
#
# Downloads dominate the first hours: 160 coverage BigWigs at about 124 MB
# each, plus one junction matrix per tissue at about 2.8 GB gzipped, all
# through the shared ETH proxy. Nothing in the workflow uses temp(), so the
# inputs stay on disk between the two configs and the second run reuses them.
#
# The working copy belongs on project storage rather than scratch: roughly
# 100 GB has to survive the weeks between the run and the figures, and scratch
# is purged after 15 days.
#
#SBATCH --job-name=fastder-gtex
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000
#SBATCH --output=slurm/logs/%x-%j.out
#SBATCH --signal=B:TERM@300

source "$(dirname "$0")/common.sh"

announce gtex-comparison gtex
echo
echo "=================== the plan ==================="
make gtex-comparison gtex "${MAKE_ARGS[@]}" EXTRA=-n 2>&1 | tail -30

echo
echo "=================== the run ==================="
# The comparison first: it is the one the paper's runtime table depends on, and
# a failure in the genome-wide atlas should not sit in front of it.
run_targets gtex-comparison
run_targets gtex
