"""Tests for data_layout.py and for the paths the rules write to.

Covers configs that share simulated reads, configs that differ in depth, seed
or chromosomes, local FASTQ input, and rules that bypass the per-config
directories.
"""

import re
from pathlib import Path

from data_layout import alignment_key, reads_key

RULES_DIR = Path(__file__).resolve().parents[1] / "workflow" / "rules"

# Directories every config shares. A rule joining one of these with a scenario
# writes where another config's scenario of the same name does.
SHARED_PREFIX = re.compile(
    r'op\.join\((DATA_DIR, "(tools|fastder)"|(LIGHT_DIR|R3_DIR), ("\{+scenario\}+"|wc\.scenario))')


def test_depth_key_ignores_number_format():
    assert reads_key(1e7, 10) == reads_key(10000000, 10) == "depth_10000000_seed_10"


def test_seeds_separate_reads():
    assert reads_key(1e7, 10) != reads_key(1e7, 11)


def test_configs_sharing_reads_and_index_share_alignments():
    base = alignment_key("config_full_simulation", "asimulator", 1e7, 10, "chr19_chr21")
    sweep = alignment_key("config_min_junction_reads_sweep", "asimulator", 1e7, 10, "chr19_chr21")
    assert base == sweep == "depth_10000000_seed_10_chr19_chr21"


def test_depth_seed_and_chromosomes_separate_alignments():
    keys = {
        alignment_key("a", "asimulator", 1e7, 10, "chr21"),
        alignment_key("a", "asimulator", 5e6, 10, "chr21"),
        alignment_key("a", "asimulator", 1e7, 11, "chr21"),
        alignment_key("a", "asimulator", 1e7, 10, "chr19_chr21"),
    }
    assert len(keys) == 4


def test_local_input_is_kept_per_config():
    assert alignment_key("config_local", "local", 0, 10, "chr21") == "config_local"


def test_no_rule_writes_a_scenario_into_a_shared_directory():
    offenders = [
        f"{path.name}:{number}"
        for path in sorted(RULES_DIR.glob("*.smk"))
        for number, line in enumerate(path.read_text().splitlines(), start=1)
        if SHARED_PREFIX.search(line)
    ]
    assert offenders == []
