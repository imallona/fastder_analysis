# Changelog

## [0.2.0] - 2026-10-16

Runs fastder 0.2.0 (`dbaf750`).

### Added

- Simulation on chr19 and chr21, ten samples, all eight ASimulatoR event classes.
- Depths of 5M, 10M, 30M and 40M reads; 10M replicates with seeds 11 and 12.
- Sweeps of `--min-junction-reads` and `--no-stitch`, and a split chain ablation.
- Coverage threshold ladders for simulation, TDP-43 and GTEx; reference 0.005 CPM.
- Run on a STAR alignment without annotation.
- fastder core scaling on genome-wide GTEx and on the simulation.
- Timed rules run one at a time, three repeats, one CPU model, host load recorded.
- Library sizes over the whole file, computed once and shared by all tools.
- Tables of tool versions, tool capabilities and every reported number.
- Supplementary figures for ablation, junction filter, scaling and GTEx transcript precision.
- Memory and runtime per rule; Slurm profile and sbatch scripts for ETH Euler.

### Fixed

- groHMM signal assigned to the wrong chromosome in multi-chromosome runs.
- `mes` labelled as mutually exclusive exons; missing event classes in plots.
- Configs overwriting each other's reads and tool outputs.
- Ablation runs included in aggregate panels.
- Disk use of uncompressed reads and kept FASTQ and BAM files.
- Reference annotation copied and prefixed once per scenario.
- Snakemake's solver loading the node during timed runs.
- Low memory requests of `figure_main_1` and the finish driver.
