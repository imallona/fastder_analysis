## fastder-evaluation entrypoints.
##
## The simulation depth sweep, the TDP-43 recount3 example, and the
## cross-depth meta-report. Each pipeline run is one config, selected through
## the FASTDER_EVAL_CONFIG environment variable (the Snakefile reads that, not
## --configfile).
##
## Usage:
##   make submodules         # fetch the fastder and monorail-external submodules
##   make sim                # the 10M paper simulation run
##   make simulations        # the full depth sweep: 5M, 10M, 30M, 40M
##   make tdp43              # TDP-43 recount3 showcase: STMN2, clean threshold
##   make tdp43-panel        # TDP-43 recount3 panel: 5 cryptic exons, low threshold
##   make gtex               # GTEx 4-tissue structural-concordance showcase
##   make gtex-smoke         # reduced GTEx run, 12 BigWigs, to validate the path
##   make meta               # render the depth-sweep report (after the runs)
##   make smoke              # quick 2-sample smoke test
##   make all                # simulations, meta, both tdp43 runs, then gtex
##   make dryrun             # snakemake -n for the 10M simulation config
##   make unlock             # release a stale snakemake lock
##
## Variables (override on the command line, e.g. make sim CORES=24):
##   CORES        snakemake --cores value (default 12)
##   ULIMIT_KB    per-process virtual memory cap in KB, inherited by every
##                job shell (default 104857600, i.e. 100 GB)
##   CONDA_ENV    conda env that holds snakemake (default snakemake)
##   CONDA_INIT   conda activation script (default ~/miniconda3/bin/activate)

CORES       ?= 12
ULIMIT_KB   ?= 104857600
CONDA_ENV   ?= snakemake
CONDA_INIT  ?= $(HOME)/miniconda3/bin/activate

WORKFLOW_DIR := workflow

## Activate the snakemake env and cap per-process virtual memory at 100 GB.
## snakemake's per-job shells inherit the ulimit, so every job is bounded.
ACTIVATE := source $(CONDA_INIT) && conda activate $(CONDA_ENV) && \
            ulimit -v $(ULIMIT_KB)

SNAKEMAKE := snakemake --cores $(CORES) -p

.DEFAULT_GOAL := help
.PHONY: help all submodules sim simulations sim-5m sim-30m sim-40m tdp43 \
        tdp43-panel gtex gtex-smoke meta smoke dryrun unlock

help:
	@echo "Targets: submodules sim simulations sim-5m sim-30m sim-40m tdp43 tdp43-panel gtex gtex-smoke meta smoke all dryrun unlock"
	@echo "Variables: CORES=$(CORES) ULIMIT_KB=$(ULIMIT_KB) CONDA_ENV=$(CONDA_ENV)"

## meta only needs the simulation results, so it runs before the tdp43 runs:
## a tdp43 failure then cannot block the cross-depth report.
all: simulations meta tdp43 tdp43-panel gtex

## Populate the git submodules. workflow/external/fastder must hold the fastder
## sources for the build_fastder rule to find a CMakeLists.txt; a plain
## git clone leaves the submodules empty. Idempotent, safe to re-run.
submodules:
	git submodule update --init --recursive

simulations: sim sim-5m sim-30m sim-40m

sim:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && \
	  FASTDER_EVAL_CONFIG=../config/config_full_simulation.yaml \
	  $(SNAKEMAKE) --use-conda --use-singularity'

sim-5m:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && \
	  FASTDER_EVAL_CONFIG=../config/config_full_simulation_5M.yaml \
	  $(SNAKEMAKE) --use-conda --use-singularity'

sim-30m:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && \
	  FASTDER_EVAL_CONFIG=../config/config_full_simulation_30M.yaml \
	  $(SNAKEMAKE) --use-conda --use-singularity'

sim-40m:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && \
	  FASTDER_EVAL_CONFIG=../config/config_full_simulation_40M.yaml \
	  $(SNAKEMAKE) --use-conda --use-singularity'

## TDP-43 recount3 showcase: a clean single threshold that isolates the STMN2
## cryptic exon. The recount3 backend has no ASimulatoR container step, so no
## --use-singularity.
tdp43:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && \
	  FASTDER_EVAL_CONFIG=../config/config_klim_2019_tdp43_recount3.yaml \
	  $(SNAKEMAKE) --use-conda'

## TDP-43 recount3 panel: a low single threshold that emits the wider cryptic
## exon panel (STMN2, HDGFL2, ELAVL3, CELF5, KCNQ2), recovered via junctions.
tdp43-panel:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && \
	  FASTDER_EVAL_CONFIG=../config/config_klim_2019_tdp43_recount3_panel.yaml \
	  $(SNAKEMAKE) --use-conda'

## GTEx structural-concordance showcase: fastder run once per tissue (brain,
## heart, skeletal muscle, whole blood) over the recount3 gtex data source,
## then the per-tissue expressed-region catalogs compared structurally. The
## recount3 backend has no ASimulatoR container step, so no --use-singularity.
gtex:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && \
	  FASTDER_EVAL_CONFIG=../config/config_gtex_concordance.yaml \
	  $(SNAKEMAKE) --use-conda'

## Reduced GTEx run: 2 tissues, 12 BigWigs, one chromosome. Exercises the
## whole gtex path cheaply; run it before make gtex to validate.
gtex-smoke:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && \
	  FASTDER_EVAL_CONFIG=../config/config_gtex_smoke.yaml \
	  $(SNAKEMAKE) --use-conda'

smoke:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && \
	  FASTDER_EVAL_CONFIG=../config/config_quick_light.yaml \
	  $(SNAKEMAKE) --use-conda --use-singularity'

## Knit the cross-depth report from the config_full_simulation* results.
## Per-run reports (summary.html, benchmarks.html, recount3.html) are produced
## by snakemake inside the run targets above.
meta:
	bash $(WORKFLOW_DIR)/scripts/render_meta_report.sh

dryrun:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && \
	  FASTDER_EVAL_CONFIG=../config/config_full_simulation.yaml \
	  snakemake --cores $(CORES) -n'

unlock:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && snakemake --unlock'
