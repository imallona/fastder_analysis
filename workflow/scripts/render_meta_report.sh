#!/usr/bin/env bash
# Render the depth-sweep meta-report from the config_full_simulation* runs.
# Reads results/config_full_simulation*/ and writes results/meta_depth_sweep.html.
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
WORKFLOW=$(dirname "$HERE")
RESULTS="$WORKFLOW/results"
OUT="$RESULTS/meta_depth_sweep.html"

rscript=""
for cand in "$WORKFLOW"/.snakemake/conda/*/bin/Rscript; do
    [ -x "$cand" ] || continue
    if "$cand" -e 'q(status = as.integer(!all(c("rmarkdown","ggplot2","dplyr","readr") %in% rownames(installed.packages()))))' 2>/dev/null; then
        rscript="$cand"
        break
    fi
done

if [ -z "$rscript" ]; then
    echo "No conda env with rmarkdown found. Run a pipeline target once so" >&2
    echo "--use-conda builds the rmarkdown environment, then rerun." >&2
    exit 1
fi

PATH="$(dirname "$rscript"):$PATH" "$rscript" -e \
    "rmarkdown::render('$WORKFLOW/reports/meta.Rmd', output_file='$OUT', \
     params = list(results_root = '$RESULTS'), quiet = TRUE)"
echo "wrote $OUT"
