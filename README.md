# FastDER Evaluation Snakemake Pipeline

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

Three quick-scope configs are included:

- `config/config.yaml`: full scope, heavy backend.
- `config/config_quick.yaml`: heavy backend, 2 samples, 100 k reads, chr21.
- `config/config_quick_light.yaml`: monorail_light backend, same 2 samples.

To pick a non-default config, set the `FASTDER_EVAL_CONFIG` environment
variable. Snakemake's own `--configfile` flag does a deep merge that unions
nested dicts like `asimulator.samples`; the env var fully replaces the
default config instead. Example:

```
FASTDER_EVAL_CONFIG=../config/config_quick_light.yaml \
  snakemake --use-conda --use-singularity --cores 12
```

### Outputs

After a successful run, `results/` contains:

- `summary.csv`: gffcompare base-level sensitivity and precision per
  (sample, parameter combination).
- `chain_stats.csv`: per-transcript statistics parsed from the fastder GTFs
  (number of exons, total exonic length, score, chromosome, strand).
- `summary.html`: rendered summary report covering the two CSVs above.
- `benchmarks.html`: rendered runtime and memory report parsed from the
  per-rule benchmark TSVs under `logs/benchmarks/`.

## Pipeline overview (not done)
`workflow/scripts/download_reference_ensembl.sh`:
This script retrieves Ensembl GRCh38 release 115 reference data, downloading chromosomes 1 through 22 and X along with the corresponding GTF, and outputs a cleaned reference set for downstream use.

## Known issues

With old conda versions, the workflow may fail during environment activation with an assertion involving `CONDA_SHLVL` / `old_conda_shlvl`. A workaround is to run `export CONDA_SHLVL=0` before starting Snakemake, or to use a recent conda version.

# Note

By default, git clone <big-repo> does not fetch the contents of submodules.

To clone including submodules: `git clone --recurse-submodules <big-repo-url>`

If you already cloned the big repo and want to populate submodules afterward: `git submodule update --init --recursive`