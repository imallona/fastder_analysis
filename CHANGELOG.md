# Changelog

## [Unreleased]

### Changed

- Library sizes come from `compute_library_sizes.py`, once per scenario. `run_derfinder.R`, `run_grohmm.R` and `run_megadepth_baseline.py` read them. Each runner computed its own before.
- Library size is now whole-file. It sums length times value, from the header. Analysed chromosomes no longer change it.
- `run_derfinder`, `run_grohmm` and `run_megadepth_baseline` declare `threads: 1`.
- `run_fastder` takes threads from `FASTDER_CORES`. `fastder.cores` sets it. `config["cores"]` is the fallback.
- Main Figure 1 drops two panels. Both plotted three tools at zero. groHMM leaves the exon accuracy and boundary panels. Its 50 nt binning cannot place exon boundaries. groHMM stays in the CDF and base-level panels. The figure runs A to F. Manuscript panel citations must move.
- Main Figure 2 drops the transcript-level panel. That level measures isoform reconstruction. No tool in the comparison attempts it. It now writes `supp_gtex_transcript_precision.pdf` instead. Nothing is deleted, so no number disappears. The figure runs A to J.
- `config_full_simulation.yaml`: `fastder.cores` is 1.
- `config_full_simulation.yaml`: chromosomes are chr21 and chr19. It was chr21 alone.
- `config_full_simulation.yaml`: ten samples per scenario. All eight ASimulatoR event classes run. It was five over four classes. Added `ir`, `a3`, `a5`, `mee` and a mixture.
- `config_gtex_comparison.yaml`: `fastder.cores` is 1.
- Corrected the `position_tolerance` comment in `config_full_simulation.yaml`. `pt0` needs exact agreement between edge and junction. It does not disable stitching.
- Depth configs are regenerated from `config_full_simulation.yaml`. All four depths share one grid.
- Combinations differing only in an ignored parameter collapse. `--no-stitch` ignores `position_tolerance` and `coverage_tolerance`.
- The parameter grid moved to `workflow/scripts/param_grid.py`. The Snakefile imports it. Job counts are unchanged.
- `run_megadepth_baseline.py` takes `--library-sizes` as optional. Without it, sizes come from the BigWig headers.

### Added

- `workflow/scripts/compute_library_sizes.py` and rule `compute_library_sizes`. Each scenario gets a `library_sizes.tsv`.
- Grid support for `min_junction_reads` (`mjr`) and `no_stitch` (`ns`). `no_stitch` is a switch. The flag is passed only when true.
- `config_full_simulation.yaml`: `no_stitch: [false, true]`.
- Rule `run_fastder_scaling`, timing fastder at each core count. It reads `fastder.scaling_cores` and `fastder.scaling_scenario`. An empty `scaling_cores` leaves it out.
- `config_gtex_concordance.yaml`: `scaling_cores: [1, 2, 4, 8, 16]`.
- `config_full_simulation.yaml`: `scaling_cores: [1, 2, 4, 8, 16]`. It contrasts with the genome-wide sweep. Parsing saturates at ten samples. Averaging saturates at two chromosomes. `make_sim_configs.py` drops the key from depth configs.
- Every executing rule declares `mem_mb` and `runtime`. Both bound concurrency locally through `--resources`. Both become scheduler requests on a cluster. The fastder rules size their request by scope. Benchmarks peaked under 2 GB per chromosome. They peaked at 30 GB genome-wide.
- Rule `record_host_info`, recording the benchmark machine. It writes CPU model, cores and memory. Output is `results/<config>/host_info.tsv`. The benchmarks report depends on it.
- `profiles/euler/config.yaml`, a profile for ETH Euler. It holds every cluster setting. The Slurm account is `es_platt`. Five timed rules pin `--constraint=EPYC_7763`. Nodes are not requested exclusively. Co-tenancy stays a caveat for Methods. Nothing under `workflow/` mentions Slurm.
- A core budget for cluster runs. Each job books its thread count. `EULER=1` passes `--resources cores_used=32`. The `es_platt` share is 208 cores. `--cores` cannot bound a cluster run. Override with `EULER_CORE_BUDGET`.
- `slurm/`, sbatch wrappers for four run groups. They are probe, simulation, GTEx and TDP-43. `slurm/README.md` covers storage and first-time setup.
- `EULER=1`, `CONDA_PREFIX_DIR`, `EXTRA` and `make envs`. `make gtex EULER=1` submits to Slurm. Without `EULER=1` nothing changes.
- `tests/test_euler_profile.py`, checking the profile against the rules. Every executing rule declares memory and runtime. The CPU pin covers exactly the timed rules. Rules using node scratch request `--tmp`.
- `scripts/collect_param_sweeps.py` and `scripts/collect_scaling.py`. They write `ablation.csv`, `min_junction_reads.csv` and `scaling.csv`. Panels plot these rather than recomputing.
- `param_grid.parse_param_id()`, the inverse of `param_id()`. A results table can be grouped by axis.
- `scripts/figures/figure_supp_revision.R` and rule `figure_supp_revision`. It writes `supp_ablation.pdf`, `supp_min_junction_reads.pdf` and `supp_scaling.pdf`. SVGs go with them. Each panel saves its data as `panel_<name>.csv`.
- The scaling panel annotates its ceilings. Parsing runs one thread per loaded sample. Averaging runs one thread per chromosome. Stitching is serial.
- `tests/test_collect_param_sweeps.py` and `tests/test_collect_scaling.py`. `test_param_grid.py` gained identifier round-trip tests. They pin the patterns `helpers.R` greps for.
- `config_min_junction_reads_sweep.yaml`, generated by `scripts/make_sim_configs.py`. It reuses the 10M simulated data. fastder runs alone, at shipped defaults. `min_junction_reads` sweeps 0, 1, 2, 5, 10, 20. Zero reproduces the published behaviour. Run it with `make mjr-sweep`.
- `scripts/figures/make_capability_table.py` and rule `capability_table`. They write `tool_capabilities.csv` and `tool_capabilities.tex`. Each cell cites its source in a comment. It replaces the two zero-bar panels.
- `tests/test_capability_table.py`, guarding the snapping cell. A bare yes there fails the test. Snapping reaches internal chain edges only.
- `tests/test_make_scenario.py`, covering the compressed round trip.

### Fixed

- Simulated reads are stored gzipped. `runASimulatoR.R` compresses them after the simulation. `make_scenario.py` writes its filtered copy compressed. STAR reads them with `--readFilesCommand zcat`. Plain FASTQ was about 680 GB. Compressed it is about 170 GB. Uncompressed reads on disk will re-simulate.
- Scenario FASTQ files and sorted BAMs are `temp()`. Nothing was reclaimed as the DAG advanced. The BAMs held 110 GB too. BigWigs, junction tables and results survive. ASimulatoR reads stay, being costly to regenerate.
- The fastder submodule tracks `revision` upstream. It tracked `main` before. The pin is `8da02f5`. It carries whole-file library size. It also carries `--min-junction-reads` and `--no-stitch`. The old pin lacked all three.
- The Makefile skips `conda activate` when `CONDA_INIT` names nothing. An already activated environment is then used as it stands. Job 11290024 died sourcing a `$HOME/miniconda3` Euler does not have.
- Defaults sit inside the checkout, so no configuration is needed. Conda environments go to `workflow/.snakemake/conda`. The container cache goes to `.apptainer-cache`. A repo on project storage keeps both off `$HOME`.
- The Euler profile sends memory as `--mem-per-cpu`. Euler's cli_filter rejects `--mem`, which `mem_mb` becomes. The plugin prefers `mem_mb_per_cpu` and emits the accepted flag. Rules keep their portable `mem_mb`; the profile translates.
- The profile pins no partition. Euler routes by requested runtime, so every rule declares one. Checked against a working run: job 11288220 pinned nothing and got `normal.24h`.
- `tests/test_euler_resources.py` keeps the translation in step. Every rule over the default declares a per-cpu figure. Totals must cover the request without wasting ten times it. Runtimes must fit their partition.
- The probe filters node features with grep. `sinfo` has no `--constraint`, so the check exited non-zero and failed the job.
- `common.sh` creates the directories it points at. apptainer refuses to build into a missing one. Job 11293212 died on the first container pull.
- `record_host_info` reads MemTotal with awk. The sed version tripped a Python SyntaxWarning on every run.
- Run paths live in `slurm/site.env`, sourced by `common.sh`. Conda environments and the image cache go to project storage. That is NFS; scratch is Lustre, which handles many small files badly. `sbatch` exports the submitting shell, so an activated environment carries over.
- The probe reads the plugin version from package metadata. Plugin 2.7.1 exposes no `__version__`.
- The CPU pin of the timed rules is the `constraint` resource. Plugin 2.7.1 sets `--constraint` itself and refuses it in `slurm_extra`. Every `run_fastder` and `run_fastder_scaling` submission of job 11304206 failed on that. `sbatch` still gets `-C EPYC_7763`.
- The sbatch wrappers populate the submodules themselves. A clone leaves them empty. `build_fastder` then has no sources. `make submodules-latest` moves a pin deliberately.
- The sbatch wrappers source `common.sh` correctly. Slurm copies the script to a spool directory. `$0` pointed there, not at the repo. Job 11285331 died after seven seconds. `SLURM_SUBMIT_DIR` locates the repo now. `common.sh` also checks its prerequisites first.
- `ml_star_align` writes scratch to `$TMPDIR`. STAR temp and sort spill move there. The Euler profile requests `--tmp` for both rules.
- Aggregate panels keep the default configuration. `--no-stitch` and `--min-junction-reads` are not accuracy settings. Folding them in moved the headline numbers. `default_grid()` in `helpers.R` holds the corner. `best_pids()` picks among default runs only.
- `meta.Rmd` filters to the default configuration. Its curves no longer average ablation runs.
- `benchmarks.Rmd` reports the machine from `host_info.tsv`. It states the sampling interval correctly. Sampling is 0.5 s, then 30 s. The text said 10 s.
- `tests/test_run_megadepth_baseline.py` used the old scoped size. It called `library_size(path, chroms)`, now gone. Rewritten for whole-file behaviour.
