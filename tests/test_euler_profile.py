"""The Euler profile can fall behind the rules it is written against.

Two things break silently. A rule added without resources runs on the profile's
4 GB default and is killed at scale. A tool added to the runtime comparison
without a CPU pin produces a wall clock measured on whatever node the scheduler
had free, which is not comparable with the others and is reported in the paper
as if it were.
"""

import re
from pathlib import Path

import pytest

# The tool envs that run this suite locally carry no PyYAML; CI installs it.
yaml = pytest.importorskip("yaml", reason="PyYAML not installed in this env")

ROOT = Path(__file__).resolve().parents[1]
RULES_DIR = ROOT / "workflow" / "rules"
PROFILE = ROOT / "profiles" / "euler" / "config.yaml"

RULE_HEADER = re.compile(r"^(?:rule|checkpoint)\s+(\w+):")

# The rules whose wall clock reaches the manuscript.
TIMED_RULES = {
    "run_fastder",
    "run_fastder_scaling",
    "run_derfinder",
    "run_grohmm",
    "run_megadepth_baseline",
}

# Aggregators: input only, no output and nothing to execute, so they never
# become a job that needs resources.
TARGET_ONLY_RULES = {"manuscript_figures"}


def rules_with_bodies():
    """Rule name to its source lines, for every rule that executes something."""
    found = {}
    for path in sorted(RULES_DIR.glob("*.smk")):
        current, body = None, []
        for line in path.read_text().splitlines():
            header = RULE_HEADER.match(line)
            if header:
                if current:
                    found[current] = body
                current, body = header.group(1), []
            elif current:
                body.append(line)
        if current:
            found[current] = body
    return {
        name: body
        for name, body in found.items()
        if name not in TARGET_ONLY_RULES
        and any(re.match(r"^    (shell|run|script|notebook|wrapper):", line) for line in body)
    }


def profile_set_resources():
    return yaml.safe_load(PROFILE.read_text())["set-resources"]


def test_every_executing_rule_declares_memory_and_runtime():
    missing = {
        name
        for name, body in rules_with_bodies().items()
        if not any("mem_mb=" in line for line in body)
        or not any("runtime=" in line for line in body)
    }
    assert not missing, sorted(missing)


def test_timed_rules_are_pinned_to_one_cpu_model():
    pinned = profile_set_resources()
    for rule in TIMED_RULES:
        extra = str(pinned.get(rule, {}).get("slurm_extra", ""))
        assert "--constraint=" in extra, rule


def test_timed_rules_do_not_take_whole_nodes():
    """--exclusive would book 128 cores per job against a 208 core share."""
    pinned = profile_set_resources()
    for rule in TIMED_RULES:
        assert "--exclusive" not in str(pinned.get(rule, {}).get("slurm_extra", "")), rule


def test_all_pinned_rules_use_the_same_model():
    """Two tools timed on different CPU models cannot be compared."""
    models = {
        model
        for resources in profile_set_resources().values()
        for model in re.findall(r"--constraint=(\S+)", str(resources.get("slurm_extra", "")))
    }
    assert len(models) == 1, models


def test_pinning_covers_exactly_the_timed_rules():
    """A tool added to the comparison without a pin fails here."""
    pinned = {
        rule
        for rule, resources in profile_set_resources().items()
        if "--constraint=" in str(resources.get("slurm_extra", ""))
    }
    assert pinned == TIMED_RULES


def test_timed_rules_exist():
    assert set(rules_with_bodies()) >= TIMED_RULES
