"""The Makefile must not insist on a conda installation that is not there.

Job 11290024 died on this: with the site paths removed from the sbatch
wrappers, nothing overrode CONDA_INIT, and make sourced a $HOME/miniconda3 that
Euler does not have while the right environment was already active.
"""

import re
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


def test_dryrun_plans_with_the_flags_of_a_run():
    result = subprocess.run(["make", "-n", "dryrun"], cwd=ROOT,
                            capture_output=True, text=True, check=True)
    assert "--use-conda" in result.stdout
    assert "--use-singularity" in result.stdout


def test_a_local_run_makes_three_passes():
    out = dry_run()
    assert out.count("snakemake --cores") == 3
    assert "--omit-from run_fastder" in out
    assert "--until run_fastder" in out
    assert "--resources timed=1" in out


def test_a_cluster_run_makes_one_pass():
    out = dry_run("EULER=1")
    assert out.count("snakemake --cores") == 1
    assert "timed=1" not in out


def test_a_dry_run_does_not_read_the_inputs():
    assert "find data/fastder" in dry_run()
    assert "find data/fastder" not in dry_run("EXTRA=-n")
    assert "find data/fastder" not in dry_run("EXTRA=--dry-run")


def listed_timed_rules():
    listed = re.search(r"^TIMED_RULES := (.+)$", (ROOT / "Makefile").read_text(), re.M)
    return listed.group(1).split()


def test_the_timed_pass_gives_every_job_the_one_slot():
    """The rules declare no such resource, so the cap alone would hold nothing back."""
    timed_pass = re.search(r"snakemake --cores \d+[^&]*--until[^&]*", dry_run()).group(0)
    assert all(f" {rule}" in timed_pass for rule in listed_timed_rules())
    assert "--default-resources timed=1" in timed_pass
    assert "--resources timed=1" in timed_pass


def test_the_timed_pass_uses_the_greedy_scheduler():
    """The default scheduler runs a solver on every core between and during jobs."""
    passes = re.findall(r"snakemake --cores \d+[^&]*", dry_run())
    assert ["--scheduler greedy" in p for p in passes] == [False, True, False]


def test_timed_rules_are_the_ones_the_profile_pins():
    """A tool added to the comparison and left out of the timed pass fails here."""
    yaml = pytest.importorskip("yaml", reason="PyYAML not installed in this env")
    profile = yaml.safe_load((ROOT / "profiles" / "euler" / "config.yaml").read_text())
    pinned = {rule for rule, resources in profile["set-resources"].items()
              if resources.get("constraint")}
    assert set(listed_timed_rules()) == pinned
