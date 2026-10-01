"""The timed rules must all repeat, or their timings are not comparable."""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RULES_DIR = ROOT / "workflow" / "rules"


def timed_rules():
    listed = re.search(r"^TIMED_RULES := (.+)$", (ROOT / "Makefile").read_text(), re.M)
    return listed.group(1).split()


def benchmark_of(rule):
    """The benchmark directive of a rule, as written."""
    for path in RULES_DIR.glob("*.smk"):
        match = re.search(rf"^rule {rule}:\n(?:    .*\n|\n)*?    benchmark:\n(.*)\n",
                          path.read_text(), re.M)
        if match:
            return match.group(1).strip()
    raise AssertionError(f"rule {rule} not found")


def test_every_timed_rule_repeats_its_benchmark():
    for rule in timed_rules():
        directive = benchmark_of(rule)
        assert directive.startswith("repeat("), rule
        assert directive.endswith("BENCHMARK_REPEATS)"), rule
