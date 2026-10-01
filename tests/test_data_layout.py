"""Tests for data_layout.py and for the paths the rules write to.

Covers configs that share simulated reads, configs that differ in depth, seed
or chromosomes, local FASTQ input, and rules that bypass the per-config
directories.
"""

import re
from pathlib import Path

from data_layout import alignment_key, index_key, reads_key, sample_seed

RULES_DIR = Path(__file__).resolve().parents[1] / "workflow" / "rules"

# Directories every config shares. A rule joining one of these with a scenario
# writes where another config's scenario of the same name does.
SHARED_PREFIX = re.compile(
    r'op\.join\((DATA_DIR, "(tools|fastder)"|(LIGHT_DIR|R3_DIR), ("\{+scenario\}+"|wc\.scenario))')


DESIGN = {"seq_depth": 10000000, "strand_specific": True,
          "samples": {"es": {"es": 1.0}, "mixed": {"es": 0.5, "ir": 0.5}}}


def test_reads_key_names_depth_and_seed():
    assert reads_key(DESIGN, 10).startswith("depth_10000000_seed_10_")


def test_depth_key_ignores_number_format():
    assert reads_key({**DESIGN, "seq_depth": 1e7}, 10) == reads_key(DESIGN, 10)


def test_seed_and_every_design_setting_separate_reads():
    keys = {
        reads_key(DESIGN, 10),
        reads_key(DESIGN, 11),
        reads_key({**DESIGN, "seq_depth": 5000000}, 10),
        reads_key({**DESIGN, "strand_specific": False}, 10),
        reads_key({**DESIGN, "samples": {"es": {"es": 1.0}}}, 10),
    }
    assert len(keys) == 5


def test_configs_sharing_reads_and_index_share_alignments():
    base = alignment_key("config_full_simulation", "asimulator", DESIGN, 10, "chr19_chr21")
    sweep = alignment_key("config_min_junction_reads_sweep", "asimulator", DESIGN, 10,
                          "chr19_chr21")
    assert base == sweep == f"{reads_key(DESIGN, 10)}_chr19_chr21"


def test_chromosomes_separate_alignments():
    assert (alignment_key("a", "asimulator", DESIGN, 10, "chr21")
            != alignment_key("a", "asimulator", DESIGN, 10, "chr19_chr21"))


def test_an_unannotated_index_separates_alignments():
    assert index_key("chr21", True) == "chr21"
    annotated = alignment_key("a", "asimulator", DESIGN, 10, index_key("chr21", True))
    unannotated = alignment_key("b", "asimulator", DESIGN, 10, index_key("chr21", False))
    assert annotated != unannotated


def test_local_input_is_kept_per_config():
    assert alignment_key("config_local", "local", None, 10, "chr21") == "config_local"


def test_every_sample_of_a_run_has_its_own_seed():
    samples = ["es", "mes", "ir"]
    seeds = [sample_seed(10, samples, sample) for sample in samples]
    assert seeds == [10000, 10001, 10002]


def test_runs_under_consecutive_seeds_share_no_sample_seed():
    samples = [f"s{i}" for i in range(10)]
    first = {sample_seed(10, samples, s) for s in samples}
    second = {sample_seed(11, samples, s) for s in samples}
    assert not first & second


def test_no_rule_writes_a_scenario_into_a_shared_directory():
    offenders = [
        f"{path.name}:{number}"
        for path in sorted(RULES_DIR.glob("*.smk"))
        for number, line in enumerate(path.read_text().splitlines(), start=1)
        if SHARED_PREFIX.search(line)
    ]
    assert offenders == []
