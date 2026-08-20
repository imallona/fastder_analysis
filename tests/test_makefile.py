"""The Makefile must not insist on a conda installation that is not there.

Job 11290024 died on this: with the site paths removed from the sbatch
wrappers, nothing overrode CONDA_INIT, and make sourced a $HOME/miniconda3 that
Euler does not have while the right environment was already active.
"""

import shutil
import subprocess

import pytest

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

pytestmark = pytest.mark.skipif(shutil.which("make") is None, reason="make not installed")


def dry_run(*overrides):
    result = subprocess.run(
        ["make", "-n", "smoke", *overrides],
        cwd=ROOT, capture_output=True, text=True, check=True,
    )
    return result.stdout


def test_a_missing_conda_init_is_skipped_rather_than_sourced():
    out = dry_run("CONDA_INIT=/nonexistent/activate")
    assert "/nonexistent/activate" not in out
    assert "conda activate" not in out
    assert "snakemake" in out


def test_an_existing_conda_init_is_still_used(tmp_path):
    fake = tmp_path / "activate"
    fake.write_text("# stand-in for a conda activate script\n")
    out = dry_run(f"CONDA_INIT={fake}", "CONDA_ENV=someenv")
    assert f"source {fake}" in out
    assert "conda activate someenv" in out


def test_euler_adds_the_profile_and_the_core_budget():
    out = dry_run("EULER=1", "CONDA_INIT=/nonexistent/activate")
    assert "--profile" in out
    assert "--resources cores_used=" in out


def test_a_local_run_adds_neither():
    out = dry_run("CONDA_INIT=/nonexistent/activate")
    assert "--profile" not in out
    assert "cores_used" not in out
