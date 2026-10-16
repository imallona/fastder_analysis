"""The sweep tables decide what the ablation and filter panels show.

The grid also moves min_coverage, min_length and position_tolerance, so a
careless average mixes runs differing for another reason. An unstitched
identifier omits position_tolerance, which an equality test would drop.
"""

import csv

import pytest

from collect_param_sweeps import (
    axis_value,
    collect,
    collect_annotation,
    comparable,
    depth_of,
    replicate_of,
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


def test_replicate_of_reads_the_suffix_and_defaults_to_one():
    assert replicate_of("/x/config_full_simulation_rep3") == 3
    assert replicate_of("/x/config_full_simulation") == 1
    assert depth_of("/x/config_full_simulation_rep3") == 10


def test_replicates_stay_separate_rows(tmp_path):
    make_run(tmp_path, "config_full_simulation",
             [summary_row("mc0.005_ml10_pt5_ns0", 60, 62)])
    make_run(tmp_path, "config_full_simulation_rep2",
             [summary_row("mc0.005_ml10_pt5_ns0", 70, 72)])
    rows = collect(str(tmp_path), "no_stitch")
    by_replicate = {r["replicate"]: r["value"] for r in rows if r["metric"] == "exon_sens"}
    assert by_replicate == {1: 60.0, 2: 70.0}
    assert {r["depth_M"] for r in rows} == {10}


def test_absent_axis_means_the_published_behaviour():
    assert axis_value({}, "min_junction_reads") == 0
    assert axis_value({}, "no_stitch") is False
    assert axis_value({"min_junction_reads": 5}, "min_junction_reads") == 5


def test_unstitched_identifiers_stay_comparable_without_position_tolerance():
    stitched = parse_param_id("mc0.005_ml10_pt5_ns0")
    unstitched = parse_param_id("mc0.005_ml10_ns1")
    assert comparable(stitched, "no_stitch")
    assert comparable(unstitched, "no_stitch")


def test_other_axes_off_default_are_excluded():
    assert not comparable(parse_param_id("mc0.01_ml10_pt5_ns0"), "no_stitch")
    assert not comparable(parse_param_id("mc0.005_ml25_pt5_ns0"), "no_stitch")
    assert not comparable(parse_param_id("mc0.005_ml10_pt20_ns0"), "no_stitch")


def test_ablation_keeps_both_arms_and_averages_over_samples(tmp_path):
    make_run(tmp_path, "config_full_simulation", [
        summary_row("mc0.005_ml10_pt5_ns0", 60, 62),
        summary_row("mc0.005_ml10_pt5_ns0", 62, 64),
        summary_row("mc0.005_ml10_ns1", 50, 52),
        # Off-default corners of the grid, which must not be folded in.
        summary_row("mc0.01_ml10_pt5_ns0", 10, 10),
        summary_row("mc0.005_ml10_pt20_ns0", 90, 90),
        summary_row("mc0.005_ml10_pt5_ns0", 99, 99, tool="derfinder"),
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
         "param_id": "mc0.005_ml10_pt5_ns0", "distance": d}
        for d in (0, 5, -5, 6, 100)
    ]
    make_run(tmp_path, "config_full_simulation",
             [summary_row("mc0.005_ml10_pt5_ns0", 60, 62)], distances)
    rows = collect(str(tmp_path), "no_stitch")
    boundary = [r for r in rows if r["metric"] == "boundary_within_5bp"][0]
    assert boundary["value"] == pytest.approx(60.0)
    assert boundary["n"] == 1


def test_junction_filter_sweep_orders_by_the_swept_value(tmp_path):
    make_run(tmp_path, "config_min_junction_reads_sweep", [
        summary_row(param_id({"min_coverage": 0.005, "min_length": 10,
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
             [summary_row("mc0.005_ml10_pt5_ns0", 60, 62)])
    rows = collect(str(tmp_path), "no_stitch")
    out = tmp_path / "ablation.csv"
    write_csv(rows, out, "no_stitch")
    with open(out) as fh:
        header = next(csv.reader(fh))
    assert header == ["depth_M", "replicate", "scenario", "tool", "no_stitch", "metric", "value", "n"]


def test_missing_results_directory_is_not_an_error(tmp_path):
    assert collect(str(tmp_path / "nothing_here"), "no_stitch") == []


def test_ablation_ignores_junction_filtered_runs(tmp_path):
    """A config sweeping both axes must not fold filtered runs into an arm."""
    make_run(tmp_path, "config_full_simulation", [
        summary_row("mc0.005_ml10_pt5_mjr0_ns0", 60, 62),
        summary_row("mc0.005_ml10_pt5_mjr20_ns0", 10, 10),
    ])
    rows = collect(str(tmp_path), "no_stitch")
    values = [r["value"] for r in rows if r["metric"] == "exon_sens"]
    assert values == [60.0]


def test_the_swept_axis_is_free_to_move(tmp_path):
    make_run(tmp_path, "config_min_junction_reads_sweep", [
        summary_row("mc0.005_ml10_pt5_mjr0_ns0", 60, 62),
        summary_row("mc0.005_ml10_pt5_mjr20_ns0", 40, 42),
    ])
    rows = collect(str(tmp_path), "min_junction_reads",
                   prefix="config_min_junction_reads_sweep")
    by_value = {r["min_junction_reads"]: r["value"]
                for r in rows if r["metric"] == "exon_sens"}
    assert by_value == {0: 60.0, 20: 40.0}


def test_annotation_table_compares_the_two_alignments_at_the_defaults(tmp_path):
    make_run(tmp_path, "config_full_simulation", [
        summary_row("mc0.005_ml10_pt5_ns0", 60, 62),
        summary_row("mc0.005_ml10_ns1", 50, 52),
        summary_row("mc0.01_ml10_pt5_ns0", 10, 10),
    ])
    make_run(tmp_path, "config_unannotated_alignment",
             [summary_row("mc0.005_ml10_pt5_ns0", 55, 61)])
    rows = collect_annotation(str(tmp_path))
    by_key = {(r["annotated_index"], r["metric"]): r["value"] for r in rows}
    assert by_key == {(1, "exon_sens"): 60.0, (1, "exon_prec"): 62.0,
                      (0, "exon_sens"): 55.0, (0, "exon_prec"): 61.0}


def test_annotation_table_has_a_row_set_per_tool(tmp_path):
    make_run(tmp_path, "config_full_simulation", [
        summary_row("mc0.005_ml10_pt5_ns0", 60, 62),
        summary_row("mc0.005_pt5", 55, 58, tool="derfinder"),
        summary_row("mc0.005_pt20", 1, 1, tool="derfinder"),
        summary_row("mc0.005", 54, 57, tool="megadepth_baseline"),
    ])
    make_run(tmp_path, "config_unannotated_alignment", [
        summary_row("mc0.005_ml10_pt5_ns0", 54, 57),
        summary_row("mc0.005_pt5", 53, 56, tool="derfinder"),
        summary_row("mc0.005", 52, 55, tool="megadepth_baseline"),
    ])
    rows = collect_annotation(str(tmp_path))
    sens = {(r["tool"], r["annotated_index"]): r["value"] for r in rows
            if r["metric"] == "exon_sens"}
    assert sens == {("fastder", 1): 60.0, ("fastder", 0): 54.0,
                    ("derfinder", 1): 55.0, ("derfinder", 0): 53.0,
                    ("megadepth_baseline", 1): 54.0, ("megadepth_baseline", 0): 52.0}


def test_the_unannotated_run_stays_out_of_the_depth_tables(tmp_path):
    make_run(tmp_path, "config_full_simulation",
             [summary_row("mc0.005_ml10_pt5_ns0", 60, 62)])
    make_run(tmp_path, "config_unannotated_alignment",
             [summary_row("mc0.005_ml10_pt5_ns0", 10, 10)])
    values = [r["value"] for r in collect(str(tmp_path), "no_stitch")
              if r["metric"] == "exon_sens"]
    assert values == [60.0]


def make_split_sample(run_dir, scenario, sample, sens, prec, distances):
    """What gffcompare and eval_fuzzy write for one split-chain sample."""
    graded = run_dir / "fastder_split" / scenario / sample / "mc0.005_ml10_pt5_ns0"
    graded.mkdir(parents=True)
    (graded / "gffcompare.stats").write_text(
        "#-----------------| Sensitivity | Precision  |\n"
        f"        Exon level:    {sens}     |    {prec}    |\n")
    write_rows(graded / "fuzzy_distances.csv", ["scenario", "sample", "param_id", "distance"],
               [{"scenario": scenario, "sample": sample,
                 "param_id": "mc0.005_ml10_pt5_ns0", "distance": d} for d in distances])


def test_ablation_has_three_configurations(tmp_path):
    run_dir = make_run(tmp_path, "config_full_simulation", [
        summary_row("mc0.005_ml10_pt5_ns0", 60, 62),
        summary_row("mc0.005_ml10_ns1", 50, 52),
    ])
    make_split_sample(run_dir, "variant_only", "es", 58.0, 61.0, (0, 3, 40, 40))
    make_split_sample(run_dir, "variant_only", "ir", 54.0, 59.0, (0, 0))
    rows = collect(str(tmp_path), "no_stitch")
    sens = {(r["tool"], r["no_stitch"]): r["value"] for r in rows if r["metric"] == "exon_sens"}
    assert sens == {("fastder", 0): 60.0, ("fastder", 1): 50.0, ("fastder_split", 0): 56.0}
    boundary = [r for r in rows
                if r["tool"] == "fastder_split" and r["metric"] == "boundary_within_5bp"]
    assert boundary[0]["value"] == pytest.approx(100.0 * 4 / 6)


def test_other_sweeps_carry_no_split_rows(tmp_path):
    run_dir = make_run(tmp_path, "config_min_junction_reads_sweep",
                       [summary_row("mc0.005_ml10_pt5_mjr0_ns0", 60, 62)])
    make_split_sample(run_dir, "variant_only", "es", 58.0, 61.0, (0,))
    rows = collect(str(tmp_path), "min_junction_reads",
                   prefix="config_min_junction_reads_sweep")
    assert {r["tool"] for r in rows} == {"fastder"}
