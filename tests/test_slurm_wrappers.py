"""The sbatch wrappers cannot use $0 to find their own directory.

Slurm copies the script to a spool directory, so $(dirname "$0") points there
and common.sh never loads. Job 11285331 died that way in seven seconds.
"""

from pathlib import Path

import pytest

SLURM_DIR = Path(__file__).resolve().parents[1] / "slurm"
WRAPPERS = sorted(p for p in SLURM_DIR.glob("*.sh") if p.name != "common.sh")


def test_there_are_wrappers_to_check():
    assert WRAPPERS


@pytest.mark.parametrize("path", WRAPPERS, ids=lambda p: p.name)
def test_wrapper_aborts_on_a_failed_source(path):
    """Without set -e a missing common.sh is a warning, not a stop."""
    text = path.read_text()
    assert "set -euo pipefail" in text


@pytest.mark.parametrize("path", WRAPPERS, ids=lambda p: p.name)
def test_wrapper_finds_the_repo_through_slurm_submit_dir(path):
    text = path.read_text()
    if "common.sh" not in text:
        pytest.skip("self-contained wrapper")
    assert "SLURM_SUBMIT_DIR" in text
    assert 'source "$(dirname "$0")/common.sh"' not in text


def test_common_sh_checks_its_prerequisites():
    text = (SLURM_DIR / "common.sh").read_text()
    for required in ("CMakeLists.txt", "command -v snakemake",
                     "snakemake_executor_plugin_slurm"):
        assert required in text


def test_scripts_populate_the_submodules_themselves():
    """A clone leaves them empty and build_fastder then has no sources."""
    for path in [SLURM_DIR / "common.sh", SLURM_DIR / "00_probe.sh"]:
        assert "git submodule update --init --recursive" in path.read_text()


def test_fastder_submodule_tracks_the_revision_branch():
    """The revision's numbers come from the revision branch of the fork."""
    modules = (SLURM_DIR.parent / ".gitmodules").read_text()
    fastder = modules.split('[submodule "workflow/external/fastder"]')[1]
    assert "url = https://github.com/imallona/fastder.git" in fastder
    assert "branch = revision" in fastder
