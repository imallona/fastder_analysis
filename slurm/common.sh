# Shared setup for the run scripts in this directory. Sourced, not executed.
#
# Every script here is one batch job that holds snakemake while snakemake
# submits each rule as its own Slurm job and waits. Nothing runs on a login
# node except environment building, which needs the proxy.

set -euo pipefail

# Slurm points TMPDIR at this node's local disk. The executor plugin submits
# child jobs with --export=ALL, so an inherited value would send them to a path
# that belongs to another node.
unset TMPDIR

# Compute nodes reach the internet only through the ETH proxy: the Ensembl
# reference, the recount3 BigWigs and junction matrices, the ASimulatoR image
# and libBigWig at cmake configure time all need it.
module load eth_proxy

CONDA_INIT="${CONDA_INIT:-/cluster/project/platt/$USER/miniforge3/bin/activate}"
CONDA_ENV="${CONDA_ENV:-snakemake}"

# Conda writes tens of thousands of small files per environment. $HOME is
# capped at 500k inodes and scratch is purged after 15 days, so the
# environments live on project storage.
CONDA_PREFIX_DIR="${CONDA_PREFIX_DIR:-/cluster/project/platt/$USER/fastder-eval-envs}"

# The container image ASimulatoR runs in, kept off $HOME for the same reason.
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-/cluster/project/platt/$USER/apptainer-cache}"

source "$CONDA_INIT"
conda activate "$CONDA_ENV"

cd "${SLURM_SUBMIT_DIR:-$PWD}"

MAKE_ARGS=(EULER=1
           CONDA_INIT="$CONDA_INIT"
           CONDA_ENV="$CONDA_ENV"
           CONDA_PREFIX_DIR="$CONDA_PREFIX_DIR")

announce() {
    echo "host $(hostname), job ${SLURM_JOB_ID:-none}, $(date -Is)"
    echo "make $* ${MAKE_ARGS[*]}"
}

# --signal reaches this shell rather than the child, so the trap forwards it
# and snakemake finishes writing its metadata instead of leaving a stale lock.
run_targets() {
    local status=0
    make "$@" "${MAKE_ARGS[@]}" &
    local child=$!
    terminate() {
        echo "wall clock approaching, stopping snakemake"
        kill -TERM "$child" 2>/dev/null || true
    }
    trap terminate TERM
    wait "$child" || status=$?
    echo "finished with status $status, $(date -Is)"
    return "$status"
}
