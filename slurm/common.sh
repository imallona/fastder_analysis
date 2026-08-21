# Shared setup for the run scripts here. Sourced, not executed.
# Each is one batch job holding snakemake while it submits and waits.

set -euo pipefail

# Child jobs inherit the environment, so this node's TMPDIR must not travel.
unset TMPDIR

# Compute nodes reach the internet only through the ETH proxy.
module load eth_proxy

cd "${SLURM_SUBMIT_DIR:-$PWD}"

# Paths for the run, kept in one file rather than spread through the scripts.
if [ -f slurm/site.env ]; then
    source slurm/site.env
fi

# sbatch exports the submitting shell, so an environment activated there is
# already on PATH. CONDA_INIT and CONDA_ENV are the fallback.
if ! command -v snakemake > /dev/null && [ -n "${CONDA_INIT:-}" ]; then
    source "$CONDA_INIT"
    conda activate "${CONDA_ENV:?CONDA_INIT is set, so CONDA_ENV must be too}"
fi

# Defaults are repo-relative, so they land wherever the checkout lives and no
# site path is needed. Snakemake's own default puts conda environments under
# workflow/.snakemake/conda, which is why CONDA_PREFIX_DIR stays empty here.
CONDA_PREFIX_DIR="${CONDA_PREFIX_DIR:-}"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-$PWD/.apptainer-cache}"
export APPTAINER_TMPDIR="${APPTAINER_TMPDIR:-$APPTAINER_CACHEDIR/tmp}"

# apptainer refuses to build into a directory that does not exist, and conda
# is no better. Job 11293212 died on the first container pull for this.
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
[ -n "$CONDA_PREFIX_DIR" ] && mkdir -p "$CONDA_PREFIX_DIR"

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
    echo "no snakemake on PATH. Activate its environment before sbatch," >&2
    echo "or set CONDA_INIT and CONDA_ENV in slurm/site.env." >&2
    exit 1
fi
if ! python -c "import snakemake_executor_plugin_slurm" 2> /dev/null; then
    echo "no Slurm executor plugin in $CONDA_ENV" >&2
    echo "  pip install snakemake-executor-plugin-slurm" >&2
    exit 1
fi

MAKE_ARGS=(EULER=1)
[ -n "${CONDA_INIT:-}" ] && MAKE_ARGS+=(CONDA_INIT="$CONDA_INIT")
[ -n "${CONDA_ENV:-}" ] && MAKE_ARGS+=(CONDA_ENV="$CONDA_ENV")
[ -n "$CONDA_PREFIX_DIR" ] && MAKE_ARGS+=(CONDA_PREFIX_DIR="$CONDA_PREFIX_DIR")

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
