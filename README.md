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

## Pipeline overview (not done)
`workflow/scripts/download_reference_ensembl.sh`:
This script retrieves Ensembl GRCh38 release 115 reference data, downloading chromosomes 1 through 22 and X along with the corresponding GTF, and outputs a cleaned reference set for downstream use.

