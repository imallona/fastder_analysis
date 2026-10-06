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


def test_a_local_run_has_a_memory_budget_in_every_pass():
    out = dry_run("CONDA_INIT=/nonexistent/activate", "MEM_MB=5000")
    assert out.count("mem_mb=5000") == 3


def test_a_cluster_run_has_no_memory_budget():
    out = dry_run("EULER=1", "CONDA_INIT=/nonexistent/activate", "MEM_MB=5000")
    assert "mem_mb=" not in out


def test_an_empty_budget_sets_no_limit():
    out = dry_run("CONDA_INIT=/nonexistent/activate", "MEM_MB=")
    assert "mem_mb=" not in out


def test_the_prepare_pass_alone_runs_no_timed_rule():
    out = dry_run("CONDA_INIT=/nonexistent/activate", "PASSES=prepare")
    assert out.count("snakemake --cores") == 1
    assert "--omit-from run_fastder" in out
    assert "--until" not in out


def test_the_timed_and_rest_passes_skip_the_preparation():
    out = dry_run("CONDA_INIT=/nonexistent/activate", "PASSES=timed rest")
    assert out.count("snakemake --cores") == 2
    assert "--omit-from" not in out
    assert "--until run_fastder" in out


def test_a_cluster_prepare_pass_omits_the_timed_rules():
    out = dry_run("EULER=1", "CONDA_INIT=/nonexistent/activate", "PASSES=prepare")
    assert out.count("snakemake --cores") == 1
    assert "--profile" in out
    assert "--omit-from run_fastder" in out


def test_euler_prepares_timed_configs_and_runs_the_others_in_full():
    result = subprocess.run(["make", "-n", "euler", "CONDA_INIT=/nonexistent/activate"],
                            cwd=ROOT, capture_output=True, text=True, check=True)
    by_config = {}
    for chunk in result.stdout.split("FASTDER_EVAL_CONFIG=../config/")[1:]:
        by_config[chunk.split(".yaml")[0]] = chunk
    assert len(by_config) == 15
    assert all("--profile" in chunk for chunk in by_config.values())
    prepared = {config for config, chunk in by_config.items() if "--omit-from run_fastder" in chunk}
    assert prepared == {
        "config_full_simulation", "config_full_simulation_5M",
        "config_full_simulation_30M", "config_full_simulation_40M",
        "config_klim_2019_tdp43_recount3", "config_klim_2019_tdp43_recount3_panel",
        "config_gtex_comparison", "config_gtex_concordance",
    }


def test_timed_configs_runs_each_timed_config_in_the_passes_asked_for():
    def recipe(*overrides):
        return subprocess.run(["make", "-n", "timed-configs", "CONDA_INIT=/nonexistent/activate", *overrides],
                              cwd=ROOT, capture_output=True, text=True, check=True).stdout
    here = recipe("PASSES=check timed", "EULER=", "QUIET_LOAD=")
    assert here.count("FASTDER_EVAL_CONFIG=../config/") == 8
    assert here.count("--until run_fastder") == 8
    assert "--profile" not in here
    assert "wait_quiet.py" not in here
    assert recipe("PASSES=rest", "EULER=1").count("--profile") == 8


def test_the_check_pass_only_plans():
    out = dry_run("CONDA_INIT=/nonexistent/activate", "PASSES=check")
    assert "check_prepared.py --commit-file data/prepared_commit.txt" in out
    assert "--root data/fastder --manifest data/prepared_manifest.tsv" in out
    assert out.count("snakemake --cores") == 1
    assert "-n || true ) | python3 scripts/check_prepared.py --plan" in out


def test_the_timed_pass_waits_for_a_quiet_machine():
    out = dry_run("CONDA_INIT=/nonexistent/activate", "QUIET_LOAD=1.5", "QUIET_WAIT_S=60")
    assert out.count("wait_quiet.py --below 1.5 --timeout 60") == 1
    assert out.index("wait_quiet.py") < out.index("--until run_fastder")


def test_an_empty_quiet_load_starts_at_once():
    assert "wait_quiet.py" not in dry_run("CONDA_INIT=/nonexistent/activate", "QUIET_LOAD=")


def test_a_snakemake_dry_run_does_not_wait():
    assert "wait_quiet.py" not in dry_run("CONDA_INIT=/nonexistent/activate", "EXTRA=-n")


def test_euler_records_the_commit_unless_it_is_a_dry_run():
    def recipe(*overrides):
        return subprocess.run(["make", "-n", "euler", "CONDA_INIT=/nonexistent/activate", *overrides],
                              cwd=ROOT, capture_output=True, text=True, check=True).stdout
    assert "git rev-parse HEAD > workflow/data/prepared_commit.txt" in recipe()
    assert "--write-manifest" in recipe()
    assert "prepared_commit.txt" not in recipe("EXTRA=-n")
    assert "--write-manifest" not in recipe("EXTRA=-n")


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


def test_the_machine_is_recorded_in_the_timed_pass():
    """Recorded earlier, its load would describe the prepare pass."""
    prepare, timed, _ = re.findall(r"snakemake --cores \d+[^&]*", dry_run())
    assert "record_host_info" in prepare.split("--omit-from")[1]
    assert "record_host_info" in timed.split("--until")[1]


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
