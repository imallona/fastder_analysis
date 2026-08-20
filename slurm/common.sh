# Shared setup for the run scripts here. Sourced, not executed.
# Each is one batch job holding snakemake while it submits and waits.

set -euo pipefail

# Child jobs inherit the environment, so this node's TMPDIR must not travel.
unset TMPDIR

# Compute nodes reach the internet only through the ETH proxy.
module load eth_proxy

CONDA_INIT="${CONDA_INIT:-/cluster/project/platt/$USER/miniforge3/bin/activate}"
CONDA_ENV="${CONDA_ENV:-snakemake}"

# $HOME caps inodes and scratch is purged, so envs live on project storage.
CONDA_PREFIX_DIR="${CONDA_PREFIX_DIR:-/cluster/project/platt/$USER/fastder-eval-envs}"

# The ASimulatoR image, off $HOME for the same reason.
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-/cluster/project/platt/$USER/apptainer-cache}"

source "$CONDA_INIT"
conda activate "$CONDA_ENV"

cd "${SLURM_SUBMIT_DIR:-$PWD}"

# build_fastder compiles these sources, and the pin is what gets benchmarked.
# A no-op when they already match.
if [ -d .git ]; then
    git submodule sync --recursive > /dev/null
    git submodule update --init --recursive
    echo "repo      $(git log --oneline -1)"
    git submodule status --recursive | sed 's/^/submodule /'
else
    echo "not a git checkout, skipping the submodule sync" >&2
fi

# Each of these cost a submission on the first Euler attempt.
if [ ! -f workflow/external/fastder/CMakeLists.txt ]; then
    echo "fastder sources are missing even after the submodule update." >&2
    echo "Check network access: compute nodes need module load eth_proxy." >&2
    exit 1
fi
if ! command -v snakemake > /dev/null; then
    echo "no snakemake after activating $CONDA_ENV from $CONDA_INIT" >&2
    echo "  conda create -n $CONDA_ENV -c conda-forge -c bioconda snakemake" >&2
    exit 1
fi
if ! python -c "import snakemake_executor_plugin_slurm" 2> /dev/null; then
    echo "no Slurm executor plugin in $CONDA_ENV" >&2
    echo "  pip install snakemake-executor-plugin-slurm" >&2
    exit 1
fi

MAKE_ARGS=(EULER=1
           CONDA_INIT="$CONDA_INIT"
           CONDA_ENV="$CONDA_ENV"
           CONDA_PREFIX_DIR="$CONDA_PREFIX_DIR")

announce() {
    echo "host $(hostname), job ${SLURM_JOB_ID:-none}, $(date -Is)"
    echo "make $* ${MAKE_ARGS[*]}"
}

# --signal hits this shell, so the trap forwards it and metadata gets written.
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
