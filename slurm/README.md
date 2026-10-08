# ETH Euler

Rules declare `mem_mb` and `runtime` only. The Slurm settings are in `profiles/euler/config.yaml` and in this directory.

## Setup

```
conda activate <snakemake env>
pip install snakemake-executor-plugin-slurm
sbatch slurm/00_probe.sh
```

- `site.env` has the paths of conda, the conda environments and the container cache, all on project storage.
- `01_prepare.sh` builds the conda environments in its own job, before it submits a rule.
- The probe checks submission from a batch job, the benchmark files and the `EPYC_7763` nodes.

## Run

```
make euler EXTRA=-n    # lists the jobs
prepare=$(sbatch --parsable slurm/01_prepare.sh)
timed=$(sbatch --parsable --dependency=afterok:$prepare slurm/02_timed.sh)
sbatch --dependency=afterok:$timed slurm/03_finish.sh
```

- `01_prepare.sh` runs the timed configs up to their timed rules and the accuracy-only configs in full, one Slurm job per rule.
- `02_timed.sh` checks the prepared tree, then runs the timed rules one at a time on its own 16 cores.
- `03_finish.sh` submits evaluation, reports and figures. It stops at a config with a timed rule left to run.
- Arguments select configs in both: `sbatch slurm/02_timed.sh tdp43 tdp43-panel`. `03_finish.sh` with arguments draws no figures.
- After a killed job: `make unlock`.

## Timing

- Timed rules: `run_fastder`, `run_fastder_scaling`, `run_derfinder`, `run_grohmm`, `run_megadepth_baseline`.
- `02_timed.sh` requests `--constraint=EPYC_7763`, the CPU model of those rules in the profile.
- Other jobs share the node. `results/<config>/host_info.tsv` has the CPU model, the cores of the allocation and the load.

## Core limit

- `EULER=1` passes `--resources cores_used=32`, and each job counts its threads. `EULER_CORE_BUDGET` sets another value.
- The profile also caps the jobs in the queue at once, `jobs: 32`. With one-core jobs this cap is reached before the core limit.
- Short jobs spend most of their time waiting in the queue, so raise both for a stage with many of them:

```
EULER_CORE_BUDGET=150 EXTRA="--jobs 150" sbatch slurm/03_finish.sh
```

- A job with more threads than this value does not start. Rules submitted by `01_prepare.sh` and `03_finish.sh` have up to 12; `run_fastder_scaling` as a Slurm job has up to 16.

## Storage

- `$HOME` has a 500k inode limit and scratch is purged after 15 days, so the clone is on project storage.
- GTEx: about 100 GB, from 160 BigWigs of 124 MB and one junction matrix per tissue.
- `SCRATCH_DIR` in `site.env`: FASTQ and BAM files go to scratch. All are `temp()`, except the 10M reads, deleted by `config_unannotated_alignment`.
- `ml_star_align` and `ml_star_index` write STAR's temporary files to `$TMPDIR`; the profile requests `--tmp` for both.
