"""Tests for choose_threshold.py.

Covers F1, the mean over samples, tools and scenarios, the choice of the best
threshold, a tie, and rows gffcompare left without a value.
"""

import csv

import pytest

from choose_threshold import choice_table, f1, ladder_table, read_ladder

COLUMNS = ["tool", "scenario", "sample", "param_id",
           "exon_sens", "exon_prec", "base_sens", "base_prec"]


def write_summary(path, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS)
        writer.writeheader()
        for tool, param_id, sample, exon_sens, exon_prec in rows:
            writer.writerow({"tool": tool, "scenario": "variant_only", "sample": sample,
                             "param_id": param_id, "exon_sens": exon_sens,
                             "exon_prec": exon_prec, "base_sens": 50, "base_prec": 50})
    return path


def test_f1_is_the_harmonic_mean():
    assert f1(50, 50) == pytest.approx(50)
    assert f1(100, 50) == pytest.approx(66.6667, rel=1e-4)
    assert f1(0, 0) == 0.0


def test_each_tool_identifier_gives_its_threshold(tmp_path):
    rows = read_ladder(write_summary(tmp_path / "s.csv", [
        ("fastder", "mc0.01_ml10_pt5_ns0", "es", 40, 60),
        ("derfinder", "mc0.01_pt5", "es", 40, 60),
        ("megadepth_baseline", "mc0.01", "es", 40, 60),
    ]))
    assert {row["min_coverage"] for row in rows} == {0.01}
    assert rows[0]["exon_f1"] == pytest.approx(48.0)


def test_the_threshold_with_the_best_mean_exon_f1_is_chosen(tmp_path):
    table = ladder_table(read_ladder(write_summary(tmp_path / "s.csv", [
        ("fastder", "mc0.005_ml10_pt5_ns0", "es", 60, 30),
        ("fastder", "mc0.005_ml10_pt5_ns0", "ir", 60, 30),
        ("derfinder", "mc0.005_pt5", "es", 60, 30),
        ("fastder", "mc0.01_ml10_pt5_ns0", "es", 50, 50),
        ("fastder", "mc0.01_ml10_pt5_ns0", "ir", 40, 60),
        ("derfinder", "mc0.01_pt5", "es", 50, 50),
        ("fastder", "mc0.05_ml10_pt5_ns0", "es", 15, 55),
        ("derfinder", "mc0.05_pt5", "es", 15, 55),
    ])))
    fastder_at_001 = [r for r in table if r["tool"] == "fastder" and r["min_coverage"] == 0.01][0]
    assert fastder_at_001["samples"] == 2
    assert fastder_at_001["exon_f1"] == pytest.approx((50 + 48) / 2)

    choice = {row["min_coverage"]: row for row in choice_table(table)}
    assert [t for t, row in choice.items() if row["chosen"]] == [0.01]
    assert choice[0.01]["tools"] == 2
    assert choice[0.01]["exon_f1"] == pytest.approx((49 + 50) / 2)


def test_a_tie_goes_to_the_higher_threshold(tmp_path):
    table = ladder_table(read_ladder(write_summary(tmp_path / "s.csv", [
        ("fastder", "mc0.01_ml10_pt5_ns0", "es", 50, 50),
        ("fastder", "mc0.02_ml10_pt5_ns0", "es", 50, 50),
    ])))
    assert [r["min_coverage"] for r in choice_table(table) if r["chosen"]] == [0.02]


def test_rows_without_a_value_are_left_out(tmp_path):
    rows = read_ladder(write_summary(tmp_path / "s.csv", [
        ("fastder", "mc0.2_ml10_pt5_ns0", "es", "", ""),
        ("grohmm", "lp-25_uts15", "es", 1, 1),
    ]))
    assert rows == []
