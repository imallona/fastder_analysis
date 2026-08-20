"""Euler's cli_filter refuses --mem, which is what mem_mb becomes.

Every job therefore needs mem_mb_per_cpu, which the plugin prefers and emits as
--mem-per-cpu. Job 11294132 died on this. The rules keep their portable mem_mb
and runtime; the profile translates, and these tests keep the two in step.
"""

import re
from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml", reason="PyYAML not installed in this env")

ROOT = Path(__file__).resolve().parents[1]
RULES_DIR = ROOT / "workflow" / "rules"
PROFILE = ROOT / "profiles" / "euler" / "config.yaml"

RULE_HEADER = re.compile(r"^(?:rule|checkpoint)\s+(\w+):")
MEM = re.compile(r"^\s*mem_mb=(\S+?),")
RUNTIME = re.compile(r"^\s*runtime=(\d+),")
THREADS = re.compile(r"^    threads:\s*(.+)$")

# The configs set cores to 12, which is what the threaded rules ask for.
CONFIGURED_CORES = 12


def declared():
    """Rule name to its declared mem_mb, runtime and threads."""
    out = {}
    for path in sorted(RULES_DIR.glob("*.smk")):
        current = None
        for line in path.read_text().splitlines():
            header = RULE_HEADER.match(line)
            if header:
                current = header.group(1)
                out[current] = {"mem_mb": None, "runtime": None, "threads": 1}
            elif current:
                mem = MEM.match(line)
                if mem:
                    out[current]["mem_mb"] = mem.group(1)
                run = RUNTIME.match(line)
                if run:
                    out[current]["runtime"] = int(run.group(1))
                threads = THREADS.match(line)
                if threads:
                    out[current]["threads"] = (
                        CONFIGURED_CORES if "cores" in threads.group(1) else threads.group(1)
                    )
    return {name: spec for name, spec in out.items() if spec["mem_mb"]}


def profile():
    return yaml.safe_load(PROFILE.read_text())


def test_defaults_cover_what_euler_requires():
    defaults = profile()["default-resources"]
    assert defaults["mem_mb_per_cpu"], "without it the plugin emits --mem"
    assert defaults["runtime"], "Euler routes to a partition by requested runtime"


def test_every_hungry_rule_has_a_per_cpu_figure():
    """A rule needing more than the default cannot rely on it."""
    default = profile()["default-resources"]["mem_mb_per_cpu"]
    overrides = profile()["set-resources"]
    for name, spec in declared().items():
        mem = spec["mem_mb"]
        if not mem.isdigit() or int(mem) <= default:
            continue
        assert "mem_mb_per_cpu" in overrides.get(name, {}), name


def test_per_cpu_figures_cover_the_declared_total():
    """mem_mb_per_cpu times threads must reach what the rule asked for."""
    overrides = profile()["set-resources"]
    for name, spec in declared().items():
        per_cpu = overrides.get(name, {}).get("mem_mb_per_cpu")
        mem = spec["mem_mb"]
        threads = spec["threads"]
        if per_cpu is None or not mem.isdigit() or not isinstance(threads, int):
            continue
        assert per_cpu * threads >= int(mem), (name, per_cpu, threads, mem)


def test_per_cpu_figures_are_not_wildly_over():
    """Ten times the request wastes a share that other people are queuing for."""
    overrides = profile()["set-resources"]
    for name, spec in declared().items():
        per_cpu = overrides.get(name, {}).get("mem_mb_per_cpu")
        mem = spec["mem_mb"]
        threads = spec["threads"]
        if per_cpu is None or not mem.isdigit() or not isinstance(threads, int):
            continue
        assert per_cpu * threads <= 10 * int(mem), (name, per_cpu, threads, mem)


def test_every_rule_declares_a_runtime():
    """Euler routes to a partition by requested runtime, so it decides where a
    job lands. platt_invivo_cortical_neurons pins no partition either, and its
    24 hour job 11288220 ran in normal.24h."""
    missing = [name for name, spec in declared().items() if spec["runtime"] is None]
    assert not missing, missing


def test_no_partition_is_pinned():
    """A pinned partition overrides that routing, and the mapping is ours to
    get wrong. The working setup on this cluster pins none."""
    overrides = profile()["set-resources"]
    assert "slurm_partition" not in profile()["default-resources"]
    assert not [r for r, v in overrides.items() if "slurm_partition" in v]


def test_the_threaded_rules_fit_one_node():
    """A 128 core node has 249 GB; asking for more queues forever."""
    overrides = profile()["set-resources"]
    for name, spec in declared().items():
        per_cpu = overrides.get(name, {}).get("mem_mb_per_cpu")
        threads = spec["threads"]
        if per_cpu is None or not isinstance(threads, int):
            continue
        assert per_cpu * threads <= 249_000, (name, per_cpu, threads)
