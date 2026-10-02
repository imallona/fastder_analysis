## Aim

Evaluation and workflow capabilities `for fastder`.

fastder enters as a git submodule ([imallona/fastder](https://github.com/imallona/fastder), a fork of [martinalavanya/fastder](https://github.com/martinalavanya/fastder)) adding lean MM/RR parsing, BigWig coverage through libBigWig, and strand-aware stitching. Clone with `git clone --recurse-submodules <url>`.

## Running

Needs Conda/Miniconda and Singularity/Apptainer. Install Snakemake into an env named `snakemake`:

```
conda create -c conda-forge -c bioconda -c nodefaults -n snakemake snakemake
make submodules        # fetch the fastder and monorail-external submodules, once after cloning
```

There is a `Makefile`; `make help` lists its targets. `make all` runs every config the figures read and then `make figures`; `make smoke` is a small end-to-end test. Override defaults on the command line, e.g. `make sim CORES=24` (`ULIMIT_KB` caps per-process virtual memory at 100 GB).

A local run makes three passes over a config: the inputs, then the tool runs whose wall clock is reported, one job at a time, then evaluation and reports. Keep the machine free of other work during the second pass.

Add `EULER=1` to any target to submit its rules to the ETH Euler cluster instead of running them here, e.g. `make gtex-comparison EULER=1`. The cluster settings live in `profiles/euler/config.yaml` and the sbatch wrappers in `slurm/`; see `slurm/README.md`. Without `EULER=1` the workflow runs locally exactly as before.

Each target runs one config; to run a config directly, set `FASTDER_EVAL_CONFIG` (it fully replaces the default config, unlike Snakemake's deep-merging `--configfile`):

```
conda activate snakemake
cd workflow/
FASTDER_EVAL_CONFIG=../config/config_full_simulation.yaml \
  snakemake --use-conda --use-singularity --cores <num_cores>
```

`--use-singularity` is required for any run with simulated input (`run_asimulator` pulls `docker://biomedbigdata/asimulator`) and for the `monorail` backend (`recount-pump`, `recount-unify`). The recount3 backend uses no container.

## Run modes / alignment backend

`monorail.backend` chooses how reads become coverage BigWigs and MM/RR junction files; everything downstream is identical.

- `monorail` (default): full Monorail stack in Singularity (STAR, BigWig, junctions, aggregation). Downloads multi-GB reference indexes on first run. Ingests fastqs.
- `monorail_light`: (perhaps chromosome-restricted) STAR, then a Python script builds lean MM/RR from `SJ.out.tab`. No Singularity, no whole-genome download. Ingests fastqs.
- `recount3`: no read processing; downloads Monorail-processed coverage and junctions recount3 already holds, reshaped per group. `recount3.data_source` is `sra` or `gtex`. Does not align.

## Configs

Our Snakemake workflow uses config files to define run properties.

- `config_full_simulation.yaml`: paper simulation, 10 samples (the eight ASimulatoR event classes and two mixtures), 10M reads, chr21 and chr19, monorail_light, 20-combination fastder grid. The 10M point of the depth sweep; `_5M`/`_30M`/`_40M` variants come from `workflow/scripts/make_sim_configs.py`.
- `config_full_simulation_rep2.yaml`, `_rep3.yaml`: the 10M simulation drawn again under seeds 11 and 12 (`make sim-replicates`). The cross-depth report and `ablation.csv` carry a `replicate` column and show the range. Written by the same script.
- `config_unannotated_alignment.yaml`: the 10M reads aligned against a STAR index built without the annotation (`monorail.annotated_index: false`), fastder alone at its defaults (`make sim-unannotated`). `annotation.csv` compares it with the annotated run. Written by the same script.
- The comparisons are made at 0.005 CPM, the threshold with the best exon-level F1 on `config_threshold_ladder.yaml`; `REFERENCE` in `workflow/scripts/param_grid.py` holds it. It is fastder's default coverage threshold; the default was 0.05 before.
- `config_threshold_ladder.yaml`: the 10M reads over eight `min_coverage` values for fastder, derfinder and the megadepth baseline (`make threshold-ladder`). `threshold_choice.csv` gives the mean exon-level and base-level F1 per threshold and flags the best. Written by the same script.
- `config_min_junction_reads_sweep.yaml`: the 10M simulated data, fastder alone at its defaults, `min_junction_reads` over 0, 1, 2, 5, 10, 20. Written by the same script.
- `config_klim_2019_tdp43_recount3.yaml`: TDP-43 knockdown vs control, motor-neuron RNA-seq (SRP166282, GSE121569), chr8/19/20. Showcase: at 0.05 CPM the STMN2 cryptic exon alone separates knockdown from control.
- `config_klim_2019_tdp43_recount3_panel.yaml`: same data at 0.005 CPM, where the most loci of the panel (STMN2, HDGFL2, ELAVL3, CELF5) separate the groups. Both thresholds come from the ladder below.
- `config_klim_2019_tdp43_recount3_ladder.yaml`: same data, fastder alone over a ladder of `min_coverage` values (`make tdp43-ladder`). `threshold_range_summary.csv` gives, per cryptic exon locus, the lowest and highest threshold at which a region is called in knockdown and not in control. The loci are in `config/tdp43_cryptic_exons.tsv`, with their source.
- `config_gtex_concordance.yaml`: fastder genome-wide on four GTEx tissues, eight sub-groups each. Clustering the 32 sub-group catalogs shows region shape carries tissue identity. `tools: [fastder]`.
- `config_gtex_threshold_ladder.yaml`: the first sub-group of each tissue on chr19, three tools over ten `min_coverage` values (`make gtex-threshold-ladder`). `threshold_choice.csv` flags the threshold with the best exon-level F1 against the annotation. Written by `workflow/scripts/make_gtex_ladder_config.py`.
- `config_gtex_comparison.yaml`: the same sub-groups on chr19 with all four tools.
- `config_local.yaml`, `config_quick(_light).yaml`, `config_medium_light.yaml`, `config.yaml`: local FASTQ and small chr21 smoke/dev runs.


### Config settings

- `fastder.chromosomes`: fastder's `--chr` and the RR filter. Omit for chr1-22 and chrX.
- `fastder.min_coverage`, `min_length`, `position_tolerance`, `coverage_tolerance`, `min_junction_reads`, `no_stitch`: lists, run as a cross-product. Omit a list for fastder's default.
- `fastder.cores`: threads for fastder itself. Defaults to `cores`; the tool comparisons set 1.
- `fastder.scaling_cores`: core counts for the fastder scaling run. Omit to skip it.
- `benchmark_repeats`: times each timed tool run is repeated. Default 1. Reports and `scaling.csv` give the median over repeats.
- `fastder.stranded`: unstranded `all.bw` vs per-strand `plus`/`minus.bw`. Not supported by the recount3 backend.
- `tools`: subset of `fastder`, `derfinder`, `megadepth_baseline`, `grohmm`. Omit to run all four.
- `asimulator.*` (when `pump_source: asimulator`): `seq_depth`, `samples` (sample to event-mix map), `probs_as_freq`, `strand_specific`.
- `seed`: the run seed. Sample number i of `asimulator.samples`, from 0, is simulated under `seed * 1000 + i`.
- `monorail.annotated_index`: `false` builds the monorail_light STAR index without `--sjdbGTFfile`. Default `true`.
- `monorail.local_samples` / `monorail.sra_samples`: for the `local` / `sra` sources.
- `recount3.data_source`, `study_acc`, `groups`: each group becomes one scenario, either a sample list under a shared `study_acc` or a `{study, samples}` map.
- `threshold_range.loci`, `case`, `control`: a loci table and two scenario names. Adds `threshold_range.csv` to the run.
- `gffcompare.reference_annotation`: truth-set annotation for real data; empty uses the downloaded reference.

## Tool comparison and params

`derfinder` (Bioconductor caller, `--cutoff`, `--min-length`, `--maxregiongap`; `workflow/scripts/run_derfinder.R`) and `megadepth_baseline` (thresholded segmenter, one transcript per run of bases at or above `--cutoff`, no stitching; `workflow/scripts/run_megadepth_baseline.py`) consume the same BigWigs. `grohmm` (HMM segmenter over 50 bp windows, `LtProbB`, `UTS`; `workflow/scripts/run_grohmm.R`) reads them too, with its own grid under `grohmm:`. Each tool writes `data/tools/<config>/{tool}/{scenario}/{param_id}/output.gtf`, graded against the same truth set (simulated GFF, or the Ensembl annotation for real data).

Shared swept parameters:

| fastder axis | megadepth_baseline | derfinder | encoded |
|---|---|---|---|
| `--min-coverage` (CPM) | `--cutoff` | `--cutoff` | `mc<v>` |
| `--min-length` (bp) | `--min-length` | `--min-length` | pinned to `fastder.min_length[0]` for baselines |
| `--position-tolerance` (bp) | (n/a) | `--maxregiongap` (analogue) | `pt<v>` (derfinder) |
| `--coverage-tolerance` | (n/a) | (n/a) | not encoded for baselines |

Grids: `fastder` is the cross-product of its config lists (`mc_ml_pt_ct_mjr_ns`), minus combinations that differ only in a parameter `--no-stitch` ignores; `derfinder` sweeps `min_coverage` x `position_tolerance` (`mc_pt`); `megadepth_baseline` sweeps `min_coverage` only (`mc`); `grohmm` sweeps `ltprobb` x `uts` (`lp_uts`). Baselines run once per (scenario, param_id) on the pooled BigWigs. To add a tool, write `run_<tool>`, add a `<tool>.yaml` env, append to `TOOLS` in `workflow/Snakefile`, and register a param-id generator in `PARAM_IDS_BY_TOOL`.
