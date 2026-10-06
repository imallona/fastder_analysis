# ETH Euler

Rules declare `mem_mb` and `runtime` only. The Slurm settings are in `profiles/euler/config.yaml` and in this directory.

## Setup

```
conda activate <snakemake env>
pip install snakemake-executor-plugin-slurm
sbatch slurm/00_probe.sh
```

- `site.env` has the paths of conda, the conda environments and the container cache, all on project storage.
- `make envs CONFIG=<config>.yaml EULER=1` builds the conda environments on a login node.
- The probe checks submission from a batch job, the benchmark files and the `EPYC_7763` nodes.

## Run

```
make euler EXTRA=-n          # lists the jobs
sbatch slurm/01_prepare.sh   # make euler
sbatch slurm/02_timed.sh     # timed rules
sbatch slurm/03_finish.sh    # evaluation, reports, figures
```

- Each script starts after the previous one has ended: `sbatch --dependency=afterok:<job id>`.
- `01_prepare.sh` runs the timed configs up to their timed rules and the accuracy-only configs in full, one Slurm job per rule.
- `02_timed.sh` checks the prepared tree, then runs the timed rules one at a time on its own 16 cores. Arguments select configs: `sbatch slurm/02_timed.sh tdp43 tdp43-panel`.
- `03_finish.sh` submits the remaining rules and the figures.
- After a killed job: `make unlock`.

## Timing

- Timed rules: `run_fastder`, `run_fastder_scaling`, `run_derfinder`, `run_grohmm`, `run_megadepth_baseline`.
- `02_timed.sh` requests `--constraint=EPYC_7763`, the CPU model of those rules in the profile.
- Other jobs share the node. `results/<config>/host_info.tsv` has the CPU model, the cores of the allocation and the load.

## Core limit

- `EULER=1` passes `--resources cores_used=32`, and each job counts its threads. `EULER_CORE_BUDGET` sets another value.
- The largest job has 12 threads and does not start under a lower value.

## Storage

- `$HOME` has a 500k inode limit and scratch is purged after 15 days, so the clone is on project storage.
- GTEx: about 100 GB, from 160 BigWigs of 124 MB and one junction matrix per tissue.
- Simulation: about 170 GB of gzipped reads. Scenario FASTQ and BAM files are `temp()`.
- `ml_star_align` and `ml_star_index` write STAR's temporary files to `$TMPDIR`; the profile requests `--tmp` for both.
