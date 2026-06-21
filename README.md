# FastDER Evaluation Snakemake Pipeline

Evaluation pipeline for fastder, accompanying the fastder paper. It runs fastder on simulated and real RNA-seq data and compares it against derfinder and a coverage-only baseline (`megadepth_baseline`). All three tools read the same coverage BigWigs and the same CPM normalisation at the same `min_coverage`, so the comparison isolates what fastder gains from splice-junction-aware stitching. Grading is `gffcompare` plus softer metrics (best Jaccard, distance to nearest boundary, locus recall, strand agreement).

Two uses:

- Method comparison on simulated data (ASimulatoR, chr21), where a known truth set gives precision and recall.
- Worked biological example: fastder on a TDP-43 knockdown and a matched control from recount3, where TDP-43 loss produces known cryptic exons.

fastder enters as a git submodule ([imallona/fastder](https://github.com/imallona/fastder), a fork of [martinalavanya/fastder](https://github.com/martinalavanya/fastder)) adding lean MM/RR parsing, BigWig coverage through libBigWig, and strand-aware stitching.

## Running

Needs Conda/Miniconda and Singularity/Apptainer. Install Snakemake into an env named `snakemake`:

```
conda create -c conda-forge -c bioconda -c nodefaults -n snakemake snakemake
make submodules        # fetch the fastder and monorail-external submodules, once after cloning
```

The root `Makefile` is the entry point; `make help` lists the targets (`sim`, `simulations`, `tdp43`, `tdp43-panel`, `gtex`, `gtex-comparison`, `gtex-pick`, `meta`, `all`, `smoke`, `dryrun`, `unlock`). Override defaults on the command line, e.g. `make sim CORES=24` (`ULIMIT_KB` caps per-process virtual memory at 100 GB).

Each target wraps one command; to run a config directly, set `FASTDER_EVAL_CONFIG` (it fully replaces the default config, unlike Snakemake's deep-merging `--configfile`):

```
conda activate snakemake
cd workflow/
FASTDER_EVAL_CONFIG=../config/config_full_simulation.yaml \
  snakemake --use-conda --use-singularity --cores <num_cores>
```

`--use-singularity` is required for any run with simulated input (`run_asimulator` pulls `docker://biomedbigdata/asimulator`) and for the `monorail` backend (`recount-pump`, `recount-unify`). The recount3 backend uses no container.

## Run modes

`monorail.backend` chooses how reads become coverage BigWigs and MM/RR junction files; everything downstream is identical.

- `monorail` (default): full Monorail stack in Singularity (STAR, BigWig, junctions, aggregation). Downloads multi-GB reference indexes on first run.
- `monorail_light`: chromosome-restricted STAR, then a Python script builds lean MM/RR from `SJ.out.tab`. No Singularity, no whole-genome download.
- `recount3`: no read processing; downloads Monorail-processed coverage and junctions recount3 already holds, reshaped per group. `recount3.data_source` is `sra` or `gtex`.

For `monorail` and `monorail_light`, `monorail.pump_source` is the read source: `asimulator`, `sra` (downloaded at run time), or `local` (paired FASTQ on disk under `monorail.local_samples`).

## Configs

- `config_full_simulation.yaml`: paper simulation, 5 samples, 10M reads, chr21, monorail_light, 8-combination fastder grid. The 10M point of the depth sweep; `_5M`/`_30M`/`_40M` variants come from `workflow/scripts/make_sim_configs.py`.
- `config_klim_2019_tdp43_recount3.yaml`: TDP-43 knockdown vs control, motor-neuron RNA-seq (SRP166282, GSE121569), chr8/19/20. Showcase: 1.0 CPM isolates the STMN2 cryptic exon.
- `config_klim_2019_tdp43_recount3_panel.yaml`: same data at 0.02 CPM so the wider panel (STMN2, HDGFL2, ELAVL3, CELF5, KCNQ2) is emitted. Only STMN2 clears the noise floor; the other four are recovered through knockdown-specific junctions. No single threshold serves both, so the example runs twice.
- `config_gtex_concordance.yaml`: fastder genome-wide on four GTEx tissues, eight sub-groups each. Clustering the 32 sub-group catalogs shows region shape carries tissue identity. `tools: [fastder]`.
- `config_gtex_comparison.yaml`: the same sub-groups on chr19 with all three tools.
- `config_local.yaml`, `config_quick(_light).yaml`, `config_medium_light.yaml`, `config.yaml`: local FASTQ and small chr21 smoke/dev runs.

## Scenarios

Each stage runs once per scenario (a `scenario` column in the results). ASimulatoR has two: `template_and_variant` (template plus alternative isoform, so a skipped exon dips) and `variant_only` (template reads dropped, so it goes to zero). recount3 input uses the sample groups (e.g. knockdown, control) as scenarios; SRA and local input use a single `all`.

## GTEx tissue set

`make gtex-pick` rewrites the `recount3.groups` block of both GTEx configs from recount3 metadata. Default `GTEX_PICK_TISSUES="BLOOD BRAIN HEART MUSCLE LIVER LUNG TESTIS ADIPOSE_TISSUE"`; other knobs `GTEX_PICK_SEED` (10), `GTEX_PICK_N` (40), `GTEX_PICK_GROUPS` (8). The helper `workflow/scripts/pick_gtex_samples.py` runs standalone.

## Outputs

`results/<config>/` holds `summary.csv` (gffcompare metrics per tool/scenario/sample/param at every level), `chain_stats.csv` (per-transcript stats from the fastder GTFs), the four `fuzzy_*.csv` (jaccard, boundary distances, locus recall, strand), and HTML reports: `summary.html`, `benchmarks.html`, plus `recount3.html` (sra source) or `gtex_concordance.html` (gtex source).

## Tool comparison

`derfinder` (Bioconductor caller, `--cutoff`, `--min-length`, `--maxregiongap`; `workflow/scripts/run_derfinder.R`) and `megadepth_baseline` (thresholded segmenter, one transcript per run of bases at or above `--cutoff`, no stitching; `workflow/scripts/run_megadepth_baseline.py`) consume the same BigWigs. Each tool writes `data/tools/{tool}/{scenario}/{param_id}/output.gtf`, graded against the same truth set (simulated GFF, or the Ensembl annotation for real data).

Shared swept parameters are encoded in `param_id` so runs are directly comparable:

| fastder axis | megadepth_baseline | derfinder | encoded |
|---|---|---|---|
| `--min-coverage` (CPM) | `--cutoff` | `--cutoff` | `mc<v>` |
| `--min-length` (bp) | `--min-length` | `--min-length` | pinned to `fastder.min_length[0]` for baselines |
| `--position-tolerance` (bp) | (n/a) | `--maxregiongap` (analogue) | `pt<v>` (derfinder) |
| `--coverage-tolerance` | (n/a) | (n/a) | not encoded for baselines |

Grids: `fastder` is the full cross-product of its four config lists (`mc_ml_pt_ct`); `derfinder` sweeps `min_coverage` x `position_tolerance` (`mc_pt`); `megadepth_baseline` sweeps `min_coverage` only (`mc`). Baselines run once per (scenario, param_id) on the pooled BigWigs. To add a tool, write `run_<tool>`, add a `<tool>.yaml` env, append to `TOOLS` in `workflow/Snakefile`, and register a param-id generator in `PARAM_IDS_BY_TOOL`.

## Config settings

- `fastder.chromosomes`: fastder's `--chr` and the RR filter. Omit for chr1-22 and chrX.
- `fastder.min_coverage`, `min_length`, `position_tolerance`, `coverage_tolerance`: lists, run as a cross-product. Omit a list for fastder's default.
- `fastder.stranded`: unstranded `all.bw` vs per-strand `plus`/`minus.bw`. Not supported by the recount3 backend.
- `tools`: subset of `fastder`, `derfinder`, `megadepth_baseline`. Omit to run all three.
- `asimulator.*` (when `pump_source: asimulator`): `seq_depth`, `samples` (sample to event-mix map), `probs_as_freq`, `strand_specific`.
- `monorail.local_samples` / `monorail.sra_samples`: for the `local` / `sra` sources.
- `recount3.data_source`, `study_acc`, `groups`: each group becomes one scenario, either a sample list under a shared `study_acc` or a `{study, samples}` map.
- `gffcompare.reference_annotation`: truth-set annotation for real data; empty uses the downloaded reference.

## Layout and tests

```
config/             per-config YAMLs
workflow/Snakefile  all rules; rules/ included files; envs/ per-rule conda yamls
workflow/scripts/   python and R helpers; reports/ Rmd templates
workflow/external/  fastder and monorail-external submodules
workflow/data/      scratch; logs/ per-rule logs and benchmarks; results/ final CSVs and HTML
tests/              pytest unit tests for workflow/scripts/ helpers
```

Run the tests from the repo root with any env that has `pytest`, `numpy`, `pyBigWig` (the `megadepth_baseline` env has all three): `pytest tests`.

## Known issues and cloning

Old conda versions can fail activation with a `CONDA_SHLVL` assertion; `export CONDA_SHLVL=0` first or use a recent conda. Clone with `git clone --recurse-submodules <url>`, or run `git submodule update --init --recursive` in an existing clone.
