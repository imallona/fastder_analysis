# Changelog

## [Unreleased]

### Changed

- Library sizes are computed once per scenario by `compute_library_sizes.py` and read by `run_derfinder.R`, `run_grohmm.R` and `run_megadepth_baseline.py`. Each runner previously computed its own, scoped to the analysed chromosomes.
- Library size is the whole-file sum of length times value, taken from the BigWig summary header. It no longer depends on which chromosomes a run analyses.
- `run_derfinder`, `run_grohmm` and `run_megadepth_baseline` declare `threads: 1`.
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

### Fixed

- `tests/test_run_megadepth_baseline.py` asserted the old chromosome-scoped library size and called `library_size(path, chroms)`, which no longer exists. Rewritten for the whole-file behaviour.
