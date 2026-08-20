## fastder-evaluation entrypoints.
##
## The simulation depth sweep, the TDP-43 recount3 example, and the
## cross-depth meta-report. Each pipeline run is one config, selected through
## the FASTDER_EVAL_CONFIG environment variable (the Snakefile reads that, not
## --configfile).
##
## Usage:
##   make submodules         # fetch the fastder and monorail-external submodules
##   make submodules-latest  # move them to the tip of the branch each tracks
##   make sim                # the 10M paper simulation run
##   make simulations        # the full depth sweep: 5M, 10M, 30M, 40M
##   make mjr-sweep          # junction read-support sensitivity, reuses the 10M simulation
##   make tdp43              # TDP-43 recount3 showcase: STMN2, clean threshold
##   make tdp43-panel        # TDP-43 recount3 panel: 5 cryptic exons, low threshold
##   make gtex               # GTEx structural-concordance atlas (genome-wide, fastder only)
##   make gtex-comparison    # GTEx tool comparison (chr19, all three tools)
##   make gtex-smoke         # reduced GTEx run, 12 BigWigs, to validate the path
##   make gtex-pick          # rewrite the gtex configs to cover the configured tissue set (BLOOD BRAIN HEART MUSCLE LIVER LUNG TESTIS ADIPOSE_TISSUE by default)
##   make meta               # render the depth-sweep report (after the runs)
##   make reports            # re-render the figure-feeding reports (PDFs included)
##   make composites         # assemble the two composite figures
##   make figures            # reports then composites
##   make smoke              # quick 2-sample smoke test
##   make all                # simulations, meta, both tdp43 runs, then gtex
##   make dryrun             # snakemake -n for the 10M simulation config
##   make unlock             # release a stale snakemake lock
##   make envs               # build every conda environment without running anything
##
## Variables also accept EXTRA, appended to every snakemake call:
##   make sim EXTRA=-n
##
## Cluster runs: add EULER=1 to any target above to submit its rules to Slurm
## through profiles/euler. slurm/ holds sbatch wrappers that do this for the
## revision's four run groups.
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

## EULER=1 sends every rule to Slurm, e.g. make gtex EULER=1.
EULER ?=
PROFILE_FLAG := $(if $(EULER),--profile $(CURDIR)/profiles/euler,)

## Cores held at once on the cluster; the es_platt share is 208.
EULER_CORE_BUDGET ?= 32
RESOURCE_FLAG := $(if $(EULER),--resources cores_used=$(EULER_CORE_BUDGET),)

## Conda env location. On Euler use project storage: $HOME caps inodes.
CONDA_PREFIX_DIR ?=
CONDA_PREFIX_FLAG := $(if $(CONDA_PREFIX_DIR),--conda-prefix $(CONDA_PREFIX_DIR),)

## Extra snakemake flags, e.g. EXTRA=-n.
EXTRA ?=

WORKFLOW_DIR := workflow

## Activate the snakemake env and cap per-process virtual memory at 100 GB.
## snakemake's per-job shells inherit the ulimit, so every job is bounded.
ACTIVATE := source $(CONDA_INIT) && conda activate $(CONDA_ENV) && \
            ulimit -v $(ULIMIT_KB)

SNAKEMAKE := snakemake --cores $(CORES) -p $(PROFILE_FLAG) $(RESOURCE_FLAG) $(CONDA_PREFIX_FLAG) $(EXTRA)

## Run snakemake targets under one config. $(1) config file, $(2) targets.
define snake
cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && FASTDER_EVAL_CONFIG=../config/$(1) $(SNAKEMAKE) --use-conda $(2)'
endef

.DEFAULT_GOAL := help
.PHONY: help all submodules sim simulations sim-5m sim-30m sim-40m tdp43 \
        tdp43-panel gtex gtex-comparison gtex-smoke gtex-pick meta reports \
        composites figures smoke dryrun unlock envs mjr-sweep submodules-latest

help:
	@echo "Targets: submodules submodules-latest sim simulations sim-5m sim-30m sim-40m mjr-sweep tdp43 tdp43-panel gtex gtex-comparison gtex-smoke gtex-pick meta reports composites figures smoke all dryrun unlock envs"
	@echo "Variables: CORES=$(CORES) ULIMIT_KB=$(ULIMIT_KB) CONDA_ENV=$(CONDA_ENV) EULER=$(EULER) EULER_CORE_BUDGET=$(EULER_CORE_BUDGET) CONDA_PREFIX_DIR=$(CONDA_PREFIX_DIR)"

## meta only needs the simulation results, so it runs before the tdp43 runs:
## a tdp43 failure then cannot block the cross-depth report.
all: simulations meta tdp43 tdp43-panel gtex gtex-comparison figures

## Submodules at their recorded commits. A clone leaves them empty.
submodules:
	git submodule update --init --recursive

## Move the pins to each tracked branch tip, fastder to revision.
submodules-latest:
	git submodule update --init --recursive --remote
	git submodule status --recursive

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

## min_junction_reads sweep, fastder alone. Reuses the 10M data, so run after sim.
mjr-sweep:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && \
	  FASTDER_EVAL_CONFIG=../config/config_min_junction_reads_sweep.yaml \
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

## GTEx structural-concordance atlas: fastder run genome-wide, once per tissue
## sub-group over the recount3 gtex data source, then the per-sub-group
## expressed-region catalogs clustered structurally. fastder only, no
## baselines. The recount3 backend has no ASimulatoR container step, so no
## --use-singularity.
gtex:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && \
	  FASTDER_EVAL_CONFIG=../config/config_gtex_concordance.yaml \
	  $(SNAKEMAKE) --use-conda'

## GTEx tool comparison: the same sub-groups on chr19 with all three tools,
## for the fastder-versus-baseline runtime and gffcompare comparison.
gtex-comparison:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && \
	  FASTDER_EVAL_CONFIG=../config/config_gtex_comparison.yaml \
	  $(SNAKEMAKE) --use-conda'

## Reduced GTEx run: 2 tissues, 12 BigWigs, one chromosome. Walks the
## whole gtex path cheaply; run it before make gtex to validate.
gtex-smoke:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && \
	  FASTDER_EVAL_CONFIG=../config/config_gtex_smoke.yaml \
	  $(SNAKEMAKE) --use-conda'

## Rewrite the recount3.groups block of both GTEx configs to cover the
## listed tissues. Existing tissues keep their sample IDs exactly (so
## already-downloaded BigWigs stay valid); missing tissues are added
## with a deterministic, seeded pick from recount3's GTEx metadata.
## Re-running is idempotent. Override TISSUES, SEED, N_PER_TISSUE or
## SUBGROUPS on the command line if needed.
GTEX_PICK_TISSUES ?= BLOOD BRAIN HEART MUSCLE LIVER LUNG TESTIS ADIPOSE_TISSUE
GTEX_PICK_SEED    ?= 10
GTEX_PICK_N       ?= 40
GTEX_PICK_GROUPS  ?= 8
GTEX_METADATA_DIR := $(WORKFLOW_DIR)/data/recount3
GTEX_RECOUNT3_URL := https://duffel.rail.bio/recount3/human/data_sources/gtex/metadata
GTEX_CONFIGS      := config/config_gtex_comparison.yaml \
                     config/config_gtex_concordance.yaml

gtex-pick:
	mkdir -p $(GTEX_METADATA_DIR)
	@for t in $(GTEX_PICK_TISSUES); do \
	  out="$(GTEX_METADATA_DIR)/$$t.recount_project.tsv"; \
	  if [ ! -f "$$out" ]; then \
	    shard=$$(printf '%s' "$$t" | tail -c 2); \
	    url="$(GTEX_RECOUNT3_URL)/$$shard/$$t/gtex.recount_project.$$t.MD.gz"; \
	    echo "fetching $$t metadata from $$url" >&2; \
	    curl -fSL "$$url" | gunzip -c > "$$out" || { echo "fetch failed for $$t" >&2; exit 1; }; \
	  fi; \
	done
	bash -c '$(ACTIVATE) && \
	  python3 $(WORKFLOW_DIR)/scripts/pick_gtex_samples.py \
	  --tissues $(GTEX_PICK_TISSUES) \
	  --metadata-dir $(GTEX_METADATA_DIR) \
	  --seed $(GTEX_PICK_SEED) \
	  --n-per-tissue $(GTEX_PICK_N) \
	  --subgroups $(GTEX_PICK_GROUPS) \
	  --apply $(GTEX_CONFIGS)'

smoke:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && \
	  FASTDER_EVAL_CONFIG=../config/config_quick_light.yaml \
	  $(SNAKEMAKE) --use-conda --use-singularity'

## Knit the cross-depth report from the config_full_simulation* results.
## Per-run reports (summary.html, benchmarks.html, recount3.html) are produced
## by snakemake inside the run targets above.
meta:
	bash $(WORKFLOW_DIR)/scripts/render_meta_report.sh

## Re-render the reports that feed the manuscript figures, emitting a vector PDF
## per figure (dev=c('png','pdf')). Cheap: rmarkdown over the existing results,
## so it needs the per-config runs already done.
reports:
	$(call snake,config_full_simulation.yaml,render_summary_report render_benchmarks_report)
	$(call snake,config_klim_2019_tdp43_recount3.yaml,render_recount3_report render_summary_report)
	$(call snake,config_gtex_concordance.yaml,render_gtex_report render_summary_report)
	$(call snake,config_gtex_comparison.yaml,render_summary_report)
	bash $(WORKFLOW_DIR)/scripts/render_meta_report.sh

## Assemble the two composite figures. The genome-wide config makes REF_GTF
## genome-wide, which the novel-exon count needs.
composites:
	$(call snake,config_gtex_concordance.yaml,manuscript_figures)

## Reports then composites.
figures: reports composites

dryrun:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && \
	  FASTDER_EVAL_CONFIG=../config/config_full_simulation.yaml \
	  snakemake --cores $(CORES) -n'

unlock:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && snakemake --unlock'

## Build one config's conda envs, on a login node: the proxy is shared.
CONFIG ?= config_full_simulation.yaml

envs:
	cd $(WORKFLOW_DIR) && bash -c '$(ACTIVATE) && \
	  FASTDER_EVAL_CONFIG=../config/$(CONFIG) \
	  $(SNAKEMAKE) --use-conda --use-singularity --conda-create-envs-only'
