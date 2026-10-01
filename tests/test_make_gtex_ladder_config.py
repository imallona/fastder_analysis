"""Tests for make_gtex_ladder_config.py.

Covers the choice of one sub-group per tissue and the settings the ladder
config changes.
"""

import pytest

yaml = pytest.importorskip("yaml", reason="PyYAML not installed in this env")

from make_gtex_ladder_config import THRESHOLD_LADDER, first_subgroups, ladder_config  # noqa: E402

BASE = {
    "cores": 12,
    "benchmark_repeats": 3,
    "recount3": {"data_source": "gtex", "groups": {
        "brain_1": {"study": "BRAIN", "samples": ["a"]},
        "brain_2": {"study": "BRAIN", "samples": ["b"]},
        "heart_1": {"study": "HEART", "samples": ["c"]},
    }},
    "fastder": {"min_coverage": [1.0], "chromosomes": ["chr19"], "cores": 1},
}


def test_one_sub_group_per_tissue_is_kept():
    assert list(first_subgroups(BASE["recount3"]["groups"])) == ["brain_1", "heart_1"]


def test_ladder_config_changes_only_what_the_ladder_needs():
    config = ladder_config(BASE)
    assert config["fastder"]["min_coverage"] == THRESHOLD_LADDER
    assert config["fastder"]["chromosomes"] == ["chr19"]
    assert config["tools"] == ["fastder", "derfinder", "megadepth_baseline"]
    assert config["threshold_choice"] is True
    assert "benchmark_repeats" not in config
    assert BASE["fastder"]["min_coverage"] == [1.0]
    assert len(BASE["recount3"]["groups"]) == 3
