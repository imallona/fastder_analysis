"""The sweep tables decide what the ablation and filter panels show.

fastder's grid also moves min_coverage, min_length and position_tolerance, so
averaging every combination into one number per swept value mixes runs that
differ for another reason. An unstitched run carries no position_tolerance in
its identifier at all, because --no-stitch makes it inert, so an equality test
on the identifier drops one arm of the ablation.
"""

import csv

import pytest

from collect_param_sweeps import (
    axis_value,
    collect,
    comparable,
    depth_of,
    write_csv,
)
from param_grid import param_id, parse_param_id

SUMMARY_COLUMNS = ["tool", "scenario", "sample", "param_id", "exon_prec", "exon_sens"]
DISTANCE_COLUMNS = ["tool", "scenario", "sample", "param_id", "distance"]


def write_rows(path, columns, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def make_run(tmp_path, name, summary_rows, distance_rows=()):
    run_dir = tmp_path / name
    run_dir.mkdir()
    write_rows(run_dir / "summary.csv", SUMMARY_COLUMNS, summary_rows)
    if distance_rows:
        write_rows(run_dir / "fuzzy_distances.csv", DISTANCE_COLUMNS, distance_rows)
    return run_dir


def summary_row(param, sens, prec, scenario="variant_only", tool="fastder"):
    return {"tool": tool, "scenario": scenario, "sample": "es", "param_id": param,
            "exon_sens": sens, "exon_prec": prec}


def test_depth_of_reads_the_suffix_and_defaults_to_ten():
    assert depth_of("/x/config_full_simulation_40M") == 40
    assert depth_of("/x/config_full_simulation") == 10


def test_absent_axis_means_the_published_behaviour():
    assert axis_value({}, "min_junction_reads") == 0
    assert axis_value({}, "no_stitch") is False
    assert axis_value({"min_junction_reads": 5}, "min_junction_reads") == 5


def test_unstitched_identifiers_stay_comparable_without_position_tolerance():
    stitched = parse_param_id("mc0.05_ml10_pt5_ns0")
    unstitched = parse_param_id("mc0.05_ml10_ns1")
    assert comparable(stitched, "no_stitch")
    assert comparable(unstitched, "no_stitch")


def test_other_axes_off_default_are_excluded():
    assert not comparable(parse_param_id("mc0.01_ml10_pt5_ns0"), "no_stitch")
    assert not comparable(parse_param_id("mc0.05_ml25_pt5_ns0"), "no_stitch")
    assert not comparable(parse_param_id("mc0.05_ml10_pt20_ns0"), "no_stitch")


def test_ablation_keeps_both_arms_and_averages_over_samples(tmp_path):
    make_run(tmp_path, "config_full_simulation", [
        summary_row("mc0.05_ml10_pt5_ns0", 60, 62),
        summary_row("mc0.05_ml10_pt5_ns0", 62, 64),
        summary_row("mc0.05_ml10_ns1", 50, 52),
        # Off-default corners of the grid, which must not be folded in.
        summary_row("mc0.01_ml10_pt5_ns0", 10, 10),
        summary_row("mc0.05_ml10_pt20_ns0", 90, 90),
        summary_row("mc0.05_ml10_pt5_ns0", 99, 99, tool="derfinder"),
    ])
    rows = collect(str(tmp_path), "no_stitch")
    by_key = {(r["no_stitch"], r["metric"]): r for r in rows}
    assert by_key[(0, "exon_sens")]["value"] == pytest.approx(61.0)
    assert by_key[(0, "exon_sens")]["n"] == 2
    assert by_key[(1, "exon_sens")]["value"] == pytest.approx(50.0)
    assert by_key[(1, "exon_prec")]["value"] == pytest.approx(52.0)


def test_boundary_share_is_a_percentage_of_boundaries(tmp_path):
    distances = [
        {"tool": "fastder", "scenario": "variant_only", "sample": "es",
         "param_id": "mc0.05_ml10_pt5_ns0", "distance": d}
        for d in (0, 5, -5, 6, 100)
    ]
    make_run(tmp_path, "config_full_simulation",
             [summary_row("mc0.05_ml10_pt5_ns0", 60, 62)], distances)
    rows = collect(str(tmp_path), "no_stitch")
    boundary = [r for r in rows if r["metric"] == "boundary_within_5bp"][0]
    assert boundary["value"] == pytest.approx(60.0)
    assert boundary["n"] == 1


def test_junction_filter_sweep_orders_by_the_swept_value(tmp_path):
    make_run(tmp_path, "config_min_junction_reads_sweep", [
        summary_row(param_id({"min_coverage": 0.05, "min_length": 10,
                              "position_tolerance": 5, "min_junction_reads": v,
                              "no_stitch": False}), 60 - v, 62 + v)
        for v in (0, 1, 5, 20)
    ])
    rows = collect(str(tmp_path), "min_junction_reads",
                   prefix="config_min_junction_reads_sweep")
    swept = sorted({r["min_junction_reads"] for r in rows})
    assert swept == [0, 1, 5, 20]


def test_written_csv_has_the_axis_as_a_column(tmp_path):
    make_run(tmp_path, "config_full_simulation",
             [summary_row("mc0.05_ml10_pt5_ns0", 60, 62)])
    rows = collect(str(tmp_path), "no_stitch")
    out = tmp_path / "ablation.csv"
    write_csv(rows, out, "no_stitch")
    with open(out) as fh:
        header = next(csv.reader(fh))
    assert header == ["depth_M", "scenario", "tool", "no_stitch", "metric", "value", "n"]


def test_missing_results_directory_is_not_an_error(tmp_path):
    assert collect(str(tmp_path / "nothing_here"), "no_stitch") == []


def test_ablation_ignores_junction_filtered_runs(tmp_path):
    """A config sweeping both axes must not fold filtered runs into an arm."""
    make_run(tmp_path, "config_full_simulation", [
        summary_row("mc0.05_ml10_pt5_mjr0_ns0", 60, 62),
        summary_row("mc0.05_ml10_pt5_mjr20_ns0", 10, 10),
    ])
    rows = collect(str(tmp_path), "no_stitch")
    values = [r["value"] for r in rows if r["metric"] == "exon_sens"]
    assert values == [60.0]


def test_the_swept_axis_is_free_to_move(tmp_path):
    make_run(tmp_path, "config_min_junction_reads_sweep", [
        summary_row("mc0.05_ml10_pt5_mjr0_ns0", 60, 62),
        summary_row("mc0.05_ml10_pt5_mjr20_ns0", 40, 42),
    ])
    rows = collect(str(tmp_path), "min_junction_reads",
                   prefix="config_min_junction_reads_sweep")
    by_value = {r["min_junction_reads"]: r["value"]
                for r in rows if r["metric"] == "exon_sens"}
    assert by_value == {0: 60.0, 20: 40.0}
