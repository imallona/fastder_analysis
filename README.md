# FastDER Evaluation Snakemake Pipeline

This pipeline simulates RNA-seq data with ASimulatoR, aligns it against chr21, and runs fastder against the simulated truth set to measure detection accuracy. The fastder code lives in a submodule pointing at [imallona/fastder](https://github.com/imallona/fastder), which is a fork of [martinalavanya/fastder](https://github.com/martinalavanya/fastder) with three extra features used by this evaluation: lean MM/RR parsing, libBigWig support so coverage is read straight from BigWig without an intermediate BedGraph, and strand-aware stitching that propagates SJ strand into the StitchedER output.

## Running the Pipeline

### Prerequisites

1. Install Conda and Snakemake. Installation instructions are available at: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html
   
   After installing Conda or Miniconda, Snakemake can be installed with:
   
   `conda create -c conda-forge -c bioconda -c nodefaults -n snakemake snakemake`

3. Install Singularity or Apptainer.

### Execution

1. Activate the Conda environment containing Snakemake:
   `conda activate snakemake`

2. Navigate to the workflow directory:
   `cd workflow/`

3. Run the pipeline:
   `snakemake --use-conda --use-singularity --cores <num_cores>`

Replace `<num_cores>` with the number of CPU cores to allocate.

`--use-singularity` is required: the `run_asimulator` rule pulls
`docker://biomedbigdata/asimulator` to provide the ASimulatoR R package, and
without that flag Snakemake skips the container directive and tries to run
the R script in the host environment, which fails with
`there is no package called 'ASimulatoR'`. The heavy backend additionally
needs Singularity for `recount-pump` and `recount-unify`.

### Backend and configs

Two backends share the same Snakefile, selected by `monorail.backend` in the
config:

- `monorail` (default): runs the full `recount-pump` and `recount-unify`
  Singularity stack. This downloads multi-GB Monorail reference indexes the
  first time it runs.
- `monorail_light`: replaces pump and unify with a chromosome-restricted STAR
  alignment plus a small Python script that emits lean MM and RR files.
  Skips the multi-hour reference download.

Five configs are included, ordered by simulation size:

- `config/config_quick_light.yaml`: 2 samples, 100k reads, chr21, monorail_light. Smoke test.
- `config/config_quick.yaml`: 2 samples, 100k reads, chr21, heavy backend.
- `config/config_medium_light.yaml`: 5 samples, 1M reads, chr21, monorail_light.
- `config/config_full_light.yaml`: 5 samples, 10M reads, chr21, monorail_light. Recommended for figures.
- `config/config.yaml`: 5 samples, 10M reads, chr21, heavy backend (full Monorail stack).

To pick a non-default config, set the `FASTDER_EVAL_CONFIG` environment
variable. Snakemake's own `--configfile` flag does a deep merge that unions
nested dicts like `asimulator.samples`; the env var fully replaces the
default config instead. Example:

```
FASTDER_EVAL_CONFIG=../config/config_quick_light.yaml \
  snakemake --use-conda --use-singularity --cores 12
```

### Scenarios

Each ASimulatoR sample is processed through two scenarios that differ in which transcripts contribute reads to the coverage track that fastder sees (a BigWig fed through libBigWig in the active build, or a BedGraph if the legacy conversion rule is wired back in):

- `template_and_variant`: ASimulatoR default. Both the canonical (template) transcript and the alternative form are simulated, so a skipped exon still has reads from the template and the coverage track shows continuous coverage with a dip rather than a hard zero.
- `variant_only`: post-filters the FASTQ to drop reads whose source transcript carries `template=TRUE` in the splicing_variants GFF, and rewrites the GFF to keep only the alternative isoforms. This isolates the AS event signal: the skipped exon has zero coverage and the gffcompare truth set contains only the alternative isoforms.

The pipeline runs every downstream stage (alignment, BigWig, MM/RR aggregation, fastder, gffcompare, fuzzy eval) once per scenario. Results carry a `scenario` column so the summary report can compare them side by side.

### Outputs

After a successful run, `results/` contains:

- `summary.csv`: gffcompare metrics per (scenario, sample, parameter combination), every level (Base, Exon, Intron, Intron chain, Transcript, Locus) plus matching counts and missed/novel ratios.
- `chain_stats.csv`: per-transcript statistics parsed from the fastder GTFs (number of exons, total exonic length, score, chromosome, strand) per (scenario, parameter combination).
- `fuzzy_jaccard.csv`: per reference transcript, the highest exonic Jaccard against any fastder transcript on the same strand. Drops gffcompare's exact-boundary requirement.
- `fuzzy_distances.csv`: signed bp distance from each fastder exon boundary to the nearest reference boundary on the same strand.
- `fuzzy_locus_recall.csv`: fraction of reference loci with at least f of their exonic length covered by any fastder exon, at thresholds f from 0.05 to 1.0 in 0.05 increments.
- `fuzzy_strand.csv`: per fastder StitchedER, classified as concordant, discordant, unstranded, or unmatched against the best-overlapping reference transcript on any strand.
- `summary.html`: rendered summary report covering all CSVs above, with a publication-style overview section faceting by AS event class.
- `benchmarks.html`: rendered runtime and memory report parsed from the per-rule benchmark TSVs under `logs/benchmarks/`.

### Reports

The rules `render_summary_report` and `render_benchmarks_report` produce `workflow/results/summary.html` and `workflow/results/benchmarks.html` at the end of a pipeline run. To rebuild only the reports without re-running upstream rules, add `--forcerun render_summary_report render_benchmarks_report` to the snakemake invocation. Both rules share the conda env at `workflow/envs/rmarkdown.yaml`, which pulls R, rmarkdown and the tidyverse plus ComplexHeatmap and circlize for the parameter-annotated heatmaps.

### Repository layout

```
config/                   per-config YAMLs, see "Backend and configs" above
workflow/Snakefile        all rules
workflow/envs/            one conda yaml per rule conda directive
workflow/scripts/         python and R helpers called by the rules
workflow/reports/         summary.Rmd and benchmarks.Rmd templates
workflow/external/        fastder and monorail-external as git submodules
data/                     pipeline scratch (asim, monorail_light, fastder)
logs/                     per-rule logs and benchmark TSVs
results/                  final CSVs and rendered HTML reports
```

### Parameters and choices

Knobs that affect what gets simulated and how fastder is evaluated:

- `asimulator.seq_depth`: total reads per sample. 10M on chr21 alone gives roughly 22x mean coverage, generous for a method-evaluation simulation.
- `asimulator.samples`: dict of sample name to AS event mix. A single key like `es: 1.0` makes every multi-exon gene carry an exon-skipping event. The `mixed` sample with four 0.25 weights distributes the four event classes across genes.
- `asimulator.probs_as_freq`: when true (our default), the values in `asimulator.samples` are treated as frequencies, so the listed values per sample sum to the fraction of multi-exon genes that carry any event. With them summing to 1.0 every gene gets one event, which is unrealistic but maximizes signal.
- `asimulator.strand_specific`: ASimulatoR is run strand-aware so the simulated reads carry strand and the truth GFF has correct strand columns.
- `monorail.backend`: `monorail` for the full recount-pump and recount-unify Singularity stack, or `monorail_light` for the lean STAR plus Python aggregator. The light backend skips the multi-hour reference download and the 2 GB of Singularity images.
- `fastder.chromosomes`: passed straight through to fastder's `--chr` flag and used to filter the RR output, so fastder only processes the named chromosomes.
- `fastder.min_coverage`, `fastder.min_length`, `fastder.position_tolerance`, `fastder.coverage_tolerance`: lists. The pipeline runs the cross-product, so each combination becomes one fastder run with its own GTF and downstream evaluation. Omit a list to use fastder's internal default.
- `fastder.stranded`: switches the BigWig pipeline between unstranded all.bw and per-strand plus.bw + minus.bw. fastder reads either via libBigWig.

## Pipeline overview (not done)
`workflow/scripts/download_reference_ensembl.sh`:
This script retrieves Ensembl GRCh38 release 115 reference data, downloading chromosomes 1 through 22 and X along with the corresponding GTF, and outputs a cleaned reference set for downstream use.

## Known issues

With old conda versions, the workflow may fail during environment activation with an assertion involving `CONDA_SHLVL` / `old_conda_shlvl`. A workaround is to run `export CONDA_SHLVL=0` before starting Snakemake, or to use a recent conda version.

# Note

By default, git clone <big-repo> does not fetch the contents of submodules.

To clone including submodules: `git clone --recurse-submodules <big-repo-url>`

If you already cloned the big repo and want to populate submodules afterward: `git submodule update --init --recursive`