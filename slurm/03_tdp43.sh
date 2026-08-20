#!/bin/bash
# The two TDP-43 runs, showcase and panel.
#
#   sbatch slurm/03_tdp43.sh
#
# Cheap next to the other two: four recount3 BigWigs over chr8, chr19 and
# chr20. It is rerun because the whole-file library size changes what a CPM
# threshold means, and the published thresholds, 1.0 CPM for the showcase and
# 0.02 CPM for the panel, were hand-tuned at the old scale.
#
# The configs still carry those old values. Re-picking them by eye a second
# time is what the revision set out to avoid, so treat this run as the input to
# the threshold sweep rather than as a finished result: it says what the
# current values now produce, and the sweep then reports the range over which
# the STMN2 cryptic exon is called in the knockdown and not in the control.
#
#SBATCH --job-name=fastder-tdp43
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2000
#SBATCH --output=slurm/logs/%x-%j.out
#SBATCH --signal=B:TERM@300

source "$(dirname "$0")/common.sh"

announce tdp43 tdp43-panel
echo
echo "=================== the plan ==================="
make tdp43 tdp43-panel "${MAKE_ARGS[@]}" EXTRA=-n 2>&1 | tail -30

echo
echo "=================== the run ==================="
run_targets tdp43
run_targets tdp43-panel
