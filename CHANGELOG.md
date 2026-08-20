# Changelog

## [Unreleased]

### Changed

- Library sizes are computed once per scenario by `compute_library_sizes.py` and read by `run_derfinder.R`, `run_grohmm.R` and `run_megadepth_baseline.py`. Each runner previously computed its own, scoped to the analysed chromosomes.
- Library size is the whole-file sum of length times value, taken from the BigWig summary header. It no longer depends on which chromosomes a run analyses.
- `run_derfinder`, `run_grohmm` and `run_megadepth_baseline` declare `threads: 1`.
- Main Figure 1 drops the multi-exon and strand-concordance panels, which plotted three tools at zero by construction, and drops groHMM from the exon-accuracy and boundary-precision panels, whose 50 nt binning cannot place an exon boundary. groHMM stays in the boundary-distance CDF and the base-level panels. The figure is now A to F rather than A to H, so panel citations in the manuscript move.
- Main Figure 2 drops the transcript-level precision panel, which measured isoform reconstruction that no tool in the comparison attempts. It is written to `supp_gtex_transcript_precision.pdf` instead of being deleted, so the number stays visible. The figure is now A to J.
- `run_fastder` takes its thread count from `FASTDER_CORES`, set by `fastder.cores` and falling back to `config["cores"]`.
- `config_full_simulation.yaml`: `fastder.cores` is 1.
- `config_full_simulation.yaml`: chromosomes are chr21 and chr19, previously chr21.
- `config_full_simulation.yaml`: ten samples per scenario over all eight ASimulatoR event classes, previously five over four classes. Added `ir`, `a3`, `a5`, `mee` and a second mixture.
- Corrected the `position_tolerance` comment in `config_full_simulation.yaml`. `pt0` requires exact coordinate agreement between an expressed region edge and a splice junction; it does not disable stitching.
- `config_gtex_comparison.yaml`: `fastder.cores` is 1.
- Depth configs regenerated from `config_full_simulation.yaml`, so all four depths share the grid.
- Parameter combinations that differ only in a parameter the run ignores are collapsed. `--no-stitch` ignores `position_tolerance` and `coverage_tolerance`.
- The parameter grid moved from the Snakefile to `workflow/scripts/param_grid.py`, which the Snakefile imports. Behaviour is unchanged: the three configs produce the same job counts as before.
- `run_megadepth_baseline.py` takes `--library-sizes` as optional. Without it the same whole-file sizes are computed from the BigWig headers, so the script still runs on its own.

### Added

- `workflow/scripts/compute_library_sizes.py` and rule `compute_library_sizes`, writing `library_sizes.tsv` per scenario.
- Grid support for `min_junction_reads` (`mjr`) and `no_stitch` (`ns`), matching the fastder flags. `no_stitch` is a switch: the flag is passed when true and omitted when false.
- `config_full_simulation.yaml`: `no_stitch: [false, true]`.
- Rule `run_fastder_scaling`, timing fastder at each core count in `fastder.scaling_cores` on `fastder.scaling_scenario`, which defaults to the first scenario. Absent or empty `scaling_cores` leaves it out of the run.
- `config_gtex_concordance.yaml`: `scaling_cores: [1, 2, 4, 8, 16]`.
- `tests/test_param_grid.py` and `tests/test_compute_library_sizes.py`.
- Every executing rule declares `mem_mb` and `runtime`. Both are plain Snakemake resources: they bound concurrency on a workstation through `--resources` and become scheduler requests on a cluster. The fastder rules size their request from the analysed scope, since the first submission's benchmarks peaked under 2 GB on one chromosome and at 30 GB genome-wide.
- Rule `record_host_info`, writing the CPU model, core count and memory of the machine that ran the benchmarks to `results/<config>/host_info.tsv`. The benchmarks report depends on it, so a wall clock is never reported without the host that produced it.
- `profiles/euler/config.yaml`, a Snakemake profile for the ETH Euler cluster. It holds every cluster specific setting, including the Slurm account and the `--constraint=EPYC_7763` pin on the five rules whose wall clock is reported in the paper. The nodes are not requested exclusively, so co-tenancy remains a caveat for Methods. Nothing under `workflow/` mentions Slurm.
- A core budget on cluster runs: every job books its thread count through `cores_used`, and `make <target> EULER=1` passes `--resources cores_used=32`, well under the 208 core `es_platt` share. `--cores` cannot do this, since Snakemake raises the local limit once jobs go to a scheduler, and `jobs` bounds jobs rather than cores. Override with `EULER_CORE_BUDGET`.
- `slurm/`, sbatch wrappers for the four run groups (probe, simulation, GTEx, TDP-43) plus `slurm/README.md` on storage layout and first-time setup.
- `EULER=1`, `CONDA_PREFIX_DIR`, `EXTRA` and a `make envs` target in the Makefile. `make gtex EULER=1` submits each rule to Slurm; without it nothing changes.
- `tests/test_euler_profile.py`, asserting that every executing rule declares memory and runtime, and that the CPU pin covers exactly the timed rules with one model across them.
- `scripts/collect_param_sweeps.py` and `scripts/collect_scaling.py`, writing `ablation.csv`, `min_junction_reads.csv` and `scaling.csv`. The panels plot these tables rather than recomputing, so a figure and the numbers on disk cannot disagree.
- `param_grid.parse_param_id()`, the inverse of `param_id()`, so a results table carrying only an identifier can be grouped by the axis that moved.
- `scripts/figures/figure_supp_revision.R` and rule `figure_supp_revision`, writing `supp_ablation.pdf`, `supp_min_junction_reads.pdf` and `supp_scaling.pdf` with their SVGs. Each panel also saves the frame it drew as `panel_<name>.csv`.
- The scaling panel annotates the structural ceilings on the wall-time facet: parsing runs one thread per loaded sample, averaging one per chromosome, and stitching is serial. An unexplained plateau reads as a limitation, a predicted one as understanding the tool.
- `tests/test_collect_param_sweeps.py`, `tests/test_collect_scaling.py`, and identifier round-trip tests in `tests/test_param_grid.py` that pin the patterns `helpers.R` greps for.
- `config_min_junction_reads_sweep.yaml`, generated by `scripts/make_sim_configs.py`: the same simulated data as the 10M run, fastder alone, every parameter at its shipped default and `min_junction_reads` swept over 0, 1, 2, 5, 10 and 20. This is the sensitivity analysis for the new junction filter; 0 reproduces the published behaviour. `make mjr-sweep`.
- `scripts/figures/make_capability_table.py` and rule `capability_table`, writing `tool_capabilities.csv` and a `tool_capabilities.tex` tabular. Every cell is read from the tool sources or the runner scripts and carries that source as a comment. It replaces the two Figure 1 panels that plotted every tool but fastder at zero.
- `tests/test_capability_table.py`, which refuses a bare "yes" in fastder's boundary-snapping cell: snapping reaches only the internal edges of a multi-exon chain, and that is the claim the reviewer objected to.
- `config_full_simulation.yaml`: `scaling_cores: [1, 2, 4, 8, 16]`, the contrasting series to the genome-wide sweep. Its ceilings are the point: parsing takes at most one thread per loaded sample and saturates at ten, averaging parallelises over chromosomes and saturates at two, and stitching is serial. `scripts/make_sim_configs.py` drops the key from the generated depth configs, since repeating one workload's sweep at four depths measures nothing new.

### Fixed

- `meta.Rmd` applies the same default-configuration filter as the figures, so the cross-depth report does not average the ablation and junction-filter runs into its curves.
- `benchmarks.Rmd` reports the machine the timings came from, reading `host_info.tsv`, and states the benchmark sampling interval correctly: every 0.5 s for the first 30 samples, then every 30 s. It previously said 10 s.
- Aggregate figure panels no longer average the ablation and junction-filter runs into the default configuration. `--no-stitch` and `--min-junction-reads` are grid axes but not accuracy settings, so folding them into the same mean moved the headline numbers with nothing in the figure saying so. `default_grid()` in `helpers.R` keeps the stitched, unfiltered corner, and `best_pids()` picks among those runs only, which stops an unstitched run winning the Jaccard and being drawn as "fastder".

- `tests/test_run_megadepth_baseline.py` asserted the old chromosome-scoped library size and called `library_size(path, chroms)`, which no longer exists. Rewritten for the whole-file behaviour.
