# Running on the ETH Euler cluster

The workflow itself knows nothing about Slurm. Rules declare `mem_mb` and `runtime`, which bound concurrency on a workstation (`snakemake --resources mem_mb=64000`) and become scheduler requests here. Everything cluster specific is in `profiles/euler/config.yaml` and in this directory; deleting both leaves the pipeline exactly as it runs locally.

`make <target> EULER=1` submits each rule as its own Slurm job instead of running it in place. The scripts here wrap that in one batch job per run group, so nothing but environment building happens on a login node.

## Where things live

`$HOME` is 50 GB with a 500k inode cap: code only. Scratch is large but purged after 15 days, which does not survive the weeks between a run and the figures it feeds. Group project storage has neither problem, so the working copy, the conda environments and the container cache belong there.

`slurm/site.env` holds those paths, in one file rather than spread through the scripts. `common.sh` sources it.

Rough sizes: the GTEx runs need about 100 GB, dominated by 160 coverage BigWigs at about 124 MB each and one junction matrix per tissue at about 2.8 GB gzipped, stored decompressed.

The simulation set is the heavy one. Measured on the first submission's run, five samples at 10M reads over chr21 produced 40 GB of FASTQ and 6.6 GB of alignments. The revision doubles the samples and runs four depths, so the same accounting gives roughly 680 GB of FASTQ and 110 GB of BAM. The reads are now stored gzipped, which takes the FASTQ side to about 170 GB.

Most of that is now transient. The scenario FASTQ files and the sorted BAMs are `temp()`: Snakemake deletes each once the jobs reading it are done, so the peak follows the concurrency rather than the total, and what survives is the coverage BigWigs, the junction tables and the results. The ASimulatoR reads themselves are kept, since regenerating them costs a full simulation. `--notemp` keeps everything, at the cost of the full footprint.

STAR's scratch and the samtools sort spill go to node-local disk through `$TMPDIR`, which is why `ml_star_align` and `ml_star_index` request `--tmp` in the profile. Only within-job temporaries can go there: anything one rule writes for another has to sit on shared storage, since the next job runs on a different node.

## First time on a new cluster

```
conda activate <the snakemake env>
pip install snakemake-executor-plugin-slurm
sbatch slurm/00_probe.sh
```

`sbatch` exports the submitting shell, so the environment activated here carries into the job. `slurm/site.env` sets where conda environments and the container image cache go; both are on project storage, which is NFS, because conda writes tens of thousands of small files per environment and Lustre handles that badly. `CONDA_INIT` and `CONDA_ENV` there are the fallback for a shell with no environment active.

Snakemake builds the per-rule conda environments in the driver before it submits anything, so there is no separate environment step. `make envs` exists for building them ahead of time, which only saves the first driver some minutes.

Environments are built on the login node on purpose. Compute nodes reach the internet only through a shared, rate limited proxy, and 32 jobs solving environments at once abuses it.

The probe settles what only a real batch job can: that the environment comes up inside a job, that a batch job may submit to Slurm, that the `benchmark:` directive still writes a TSV when the payload runs on a compute node, and that the `EPYC_7763` nodes the profile pins to are reachable from this account.

## The runs

```
sbatch slurm/01_simulation.sh    # four depths, eight event classes, chr21 and chr19, the --no-stitch ablation
sbatch slurm/02_gtex.sh          # chr19 four-tool comparison, then the genome-wide atlas and the scaling sweep
sbatch slurm/03_tdp43.sh         # showcase and panel
```

Each is one driver job holding snakemake. `--signal=B:TERM@300` gives it five minutes to stop cleanly before the wall clock, so it writes its metadata instead of leaving a stale lock. If a driver is killed outright, `make unlock` clears the lock.

## Why the timed rules are pinned

`run_fastder`, `run_fastder_scaling`, `run_derfinder`, `run_grohmm` and `run_megadepth_baseline` are the rules whose wall clock is reported in the paper. `profiles/euler/config.yaml` gives them `--constraint=EPYC_7763`.

Euler's normal partitions mix EPYC 9654, 7742, 7H12 and 7763, and a single-core wall clock measured on one generation is not comparable with the same run on another. EPYC_7763 has 246 nodes, so availability is good.

The nodes are not requested exclusively. An exclusive job holds all 128 cores whatever its thread count, so 32 concurrent timed jobs would book 4096 cores against a 208 core share. The cost is co-tenancy: fastder is memory bandwidth bound while parsing, so a neighbour on the other sockets can inflate a wall clock. That belongs in Methods as a caveat.

`make <target> EULER=1` also passes `--resources cores_used=32`, well under the 208 core `es_platt` share, so this benchmark does not crowd out the rest of the group. Each job books its thread count, which is what bounds the workflow: `--cores` does not, since Snakemake raises the local limit once jobs go to a scheduler, and `jobs` counts jobs rather than cores. At 32 cores that is 32 single-core tool runs at a time, or two of the twelve-thread rules. Raise it with `EULER_CORE_BUDGET` if a run needs to finish sooner.

The budget has to stay at or above the largest single job. The widest are the twelve-thread rules and the top point of the scaling sweep, `fastder.scaling_cores`, which is 16. A job asking for more cores than the budget never becomes runnable.

`record_host_info` writes the CPU model, core count and memory of the machine that ran the benchmarks into `results/<config>/host_info.tsv`, and the benchmarks report depends on it. That is where the Methods sentence about the benchmark machine comes from, rather than from memory.

`tests/test_euler_profile.py` asserts that the pinned set still matches the timed rules, so adding a tool to the comparison without pinning it fails the test rather than quietly producing an incomparable number.
