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

### Added

- `workflow/scripts/compute_library_sizes.py` and rule `compute_library_sizes`, writing `library_sizes.tsv` per scenario.
