## Aim

Evaluation and workflow capabilities `for fastder`.

fastder enters as a git submodule ([imallona/fastder](https://github.com/imallona/fastder), a fork of [martinalavanya/fastder](https://github.com/martinalavanya/fastder)) adding lean MM/RR parsing, BigWig coverage through libBigWig, and strand-aware stitching. Clone with `git clone --recurse-submodules <url>`.

## Running

Needs Conda/Miniconda and Singularity/Apptainer. Install Snakemake into an env named `snakemake`:

```
conda create -c conda-forge -c bioconda -c nodefaults -n snakemake snakemake
make submodules        # fetch the fastder and monorail-external submodules, once after cloning
```

There is a `Makefile` including `make help` (`sim`, `simulations`, `tdp43`, `tdp43-panel`, `gtex`, `gtex-comparison`, `gtex-pick`, `meta`, `all`, `smoke`, `dryrun`, `unlock`). Override defaults on the command line, e.g. `make sim CORES=24` (`ULIMIT_KB` caps per-process virtual memory at 100 GB).

Each target wraps one command; to run a config directly, set `FASTDER_EVAL_CONFIG` (it fully replaces the default config, unlike Snakemake's deep-merging `--configfile`):

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

- `config_full_simulation.yaml`: paper simulation, 5 samples, 10M reads, chr21, monorail_light, 8-combination fastder grid. The 10M point of the depth sweep; `_5M`/`_30M`/`_40M` variants come from `workflow/scripts/make_sim_configs.py`.
- `config_klim_2019_tdp43_recount3.yaml`: TDP-43 knockdown vs control, motor-neuron RNA-seq (SRP166282, GSE121569), chr8/19/20. Showcase: 1.0 CPM isolates the STMN2 cryptic exon.
- `config_klim_2019_tdp43_recount3_panel.yaml`: same data at 0.02 CPM so the wider panel (STMN2, HDGFL2, ELAVL3, CELF5, KCNQ2) is emitted. Only STMN2 clears the noise floor; the other four are recovered through knockdown-specific junctions. No single threshold serves both, so the example runs twice.
- `config_gtex_concordance.yaml`: fastder genome-wide on four GTEx tissues, eight sub-groups each. Clustering the 32 sub-group catalogs shows region shape carries tissue identity. `tools: [fastder]`.
- `config_gtex_comparison.yaml`: the same sub-groups on chr19 with all three tools.
- `config_local.yaml`, `config_quick(_light).yaml`, `config_medium_light.yaml`, `config.yaml`: local FASTQ and small chr21 smoke/dev runs.


### Config settings

- `fastder.chromosomes`: fastder's `--chr` and the RR filter. Omit for chr1-22 and chrX.
- `fastder.min_coverage`, `min_length`, `position_tolerance`, `coverage_tolerance`: lists, run as a cross-product. Omit a list for fastder's default.
- `fastder.stranded`: unstranded `all.bw` vs per-strand `plus`/`minus.bw`. Not supported by the recount3 backend.
- `tools`: subset of `fastder`, `derfinder`, `megadepth_baseline`. Omit to run all three.
- `asimulator.*` (when `pump_source: asimulator`): `seq_depth`, `samples` (sample to event-mix map), `probs_as_freq`, `strand_specific`.
- `monorail.local_samples` / `monorail.sra_samples`: for the `local` / `sra` sources.
- `recount3.data_source`, `study_acc`, `groups`: each group becomes one scenario, either a sample list under a shared `study_acc` or a `{study, samples}` map.
- `gffcompare.reference_annotation`: truth-set annotation for real data; empty uses the downloaded reference.

## Tool comparison and params

`derfinder` (Bioconductor caller, `--cutoff`, `--min-length`, `--maxregiongap`; `workflow/scripts/run_derfinder.R`) and `megadepth_baseline` (thresholded segmenter, one transcript per run of bases at or above `--cutoff`, no stitching; `workflow/scripts/run_megadepth_baseline.py`) consume the same BigWigs. Each tool writes `data/tools/{tool}/{scenario}/{param_id}/output.gtf`, graded against the same truth set (simulated GFF, or the Ensembl annotation for real data).

Shared swept parameters:

| fastder axis | megadepth_baseline | derfinder | encoded |
|---|---|---|---|
| `--min-coverage` (CPM) | `--cutoff` | `--cutoff` | `mc<v>` |
| `--min-length` (bp) | `--min-length` | `--min-length` | pinned to `fastder.min_length[0]` for baselines |
| `--position-tolerance` (bp) | (n/a) | `--maxregiongap` (analogue) | `pt<v>` (derfinder) |
| `--coverage-tolerance` | (n/a) | (n/a) | not encoded for baselines |

Grids: `fastder` is the full cross-product of its four config lists (`mc_ml_pt_ct`); `derfinder` sweeps `min_coverage` x `position_tolerance` (`mc_pt`); `megadepth_baseline` sweeps `min_coverage` only (`mc`). Baselines run once per (scenario, param_id) on the pooled BigWigs. To add a tool, write `run_<tool>`, add a `<tool>.yaml` env, append to `TOOLS` in `workflow/Snakefile`, and register a param-id generator in `PARAM_IDS_BY_TOOL`.
