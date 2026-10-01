"""Tests for collect_reported_numbers.py.

Builds a small results and benchmark tree and checks the default parameters
picked per tool, each kind of value, repeats and replicates, a missing input,
and the two written files.
"""

import csv

import pytest

from collect_reported_numbers import (
    benchmark_medians,
    collect,
    reference_params,
    run_label,
    tex_value,
    write_csv,
    write_tex,
)

FASTDER = "mc0.005_ml10_pt5_ns0"
SCENARIO = "variant_only"


def write_rows(path, rows, delimiter=","):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def make_run(results_root, config="config_full_simulation"):
    run_dir = results_root / config

    def row(tool, param_id, sample="es", **values):
        return {"tool": tool, "scenario": SCENARIO, "sample": sample, "param_id": param_id,
                **values}

    write_rows(run_dir / "summary.csv", [
        row("fastder", FASTDER, "es", exon_prec=60, exon_sens=40, base_prec=90, base_sens=""),
        row("fastder", FASTDER, "ir", exon_prec=50, exon_sens=30, base_prec=92, base_sens=""),
        row("fastder", "mc0.005_ml10_ns1", "es", exon_prec=1, exon_sens=1, base_prec=1, base_sens=""),
        row("fastder", "mc0.01_ml10_pt5_ns0", "es", exon_prec=2, exon_sens=2, base_prec=2,
            base_sens=""),
        row("derfinder", "mc0.005_pt5", "es", exon_prec=55, exon_sens=41, base_prec=96,
            base_sens=""),
        row("derfinder", "mc0.005_pt20", "es", exon_prec=3, exon_sens=3, base_prec=3, base_sens=""),
        row("grohmm", "lp-25_uts15", "es", exon_prec=0.2, exon_sens=0.1, base_prec=74,
            base_sens=""),
        row("grohmm", "lp-200_uts10", "es", exon_prec=9, exon_sens=9, base_prec=9, base_sens=""),
    ])
    write_rows(run_dir / "fuzzy_jaccard.csv", [
        row("fastder", FASTDER, jaccard=0.2), row("fastder", FASTDER, jaccard=0.4),
        row("fastder", FASTDER, jaccard=0.9), row("fastder", "mc0.005_ml10_ns1", jaccard=1.0),
        row("grohmm", "lp-25_uts15", jaccard=0.5), row("grohmm", "lp-200_uts10", jaccard=0.1),
    ])
    write_rows(run_dir / "fuzzy_distances.csv", [
        row("fastder", FASTDER, "es", distance=d) for d in (0, 5, -40, -50)
    ] + [row("fastder", FASTDER, "ir", distance=0)])
    write_rows(run_dir / "fuzzy_locus_recall.csv", [
        row("fastder", FASTDER, "es", threshold=0.5, recall=0.4),
        row("fastder", FASTDER, "ir", threshold=0.5, recall=0.6),
        row("fastder", FASTDER, "es", threshold=0.9, recall=0.1),
    ])
    write_rows(run_dir / "fuzzy_strand.csv", [
        row("fastder", FASTDER, category=category, n_fastder_transcripts=count)
        for category, count in (("concordant", 30), ("discordant", 10),
                                ("unstranded", 55), ("unmatched", 5))
    ])
    return run_dir


def make_benchmark(bench_root, config, rule, job, repeats):
    """A benchmark file with one row per repeat, each a (wall, rss) pair."""
    write_rows(bench_root / config / rule / f"{job}.tsv",
               [{"s": wall, "max_rss": rss} for wall, rss in repeats], delimiter="\t")


@pytest.fixture
def tree(tmp_path):
    results, bench = tmp_path / "results", tmp_path / "bench"
    make_run(results)
    config = "config_full_simulation"
    make_benchmark(bench, config, "run_fastder", "a", [(10, 100), (12, 110), (50, 120)])
    make_benchmark(bench, config, "run_fastder", "b", [(8, ""), (8, ""), (8, "")])
    make_benchmark(bench, config, "run_derfinder", "a", [(100, 4000), (100, 4000), (100, 4000)])
    make_benchmark(bench, config, "run_fastder_scaling", "cores4", [(6, 2048)])
    return results, bench


def values(tree, **kwargs):
    results, bench = tree
    numbers = collect(str(results), str(bench), ["config_full_simulation"], **kwargs)
    return {n["name"]: n["value"] for n in numbers}


def test_run_label_names_depth_and_replicate():
    assert run_label("config_full_simulation") == "sim.10M"
    assert run_label("config_full_simulation_40M") == "sim.40M"
    assert run_label("config_full_simulation_rep3") == "sim.10M.rep3"


def test_default_parameters_per_tool(tree):
    params = reference_params(str(tree[0] / "config_full_simulation"))
    assert params == {"fastder": FASTDER, "derfinder": "mc0.005_pt5",
                      "megadepth_baseline": "mc0.005", "grohmm": "lp-25_uts15"}


def test_accuracy_is_the_mean_over_samples_at_the_defaults(tree):
    found = values(tree)
    assert found[f"sim.10M.{SCENARIO}.fastder.exon_prec"] == pytest.approx(55.0)
    assert found[f"sim.10M.{SCENARIO}.derfinder.exon_prec"] == pytest.approx(55.0)
    assert found[f"sim.10M.{SCENARIO}.grohmm.exon_prec"] == pytest.approx(0.2)
    assert found[f"sim.10M.{SCENARIO}.fastder.ir.exon_sens"] == pytest.approx(30.0)
    assert f"sim.10M.{SCENARIO}.fastder.base_sens" not in found


def test_fuzzy_values(tree):
    found = values(tree)
    assert found[f"sim.10M.{SCENARIO}.fastder.jaccard_median"] == pytest.approx(0.4)
    # es has two of four boundaries within 5 bp, ir one of one.
    assert found[f"sim.10M.{SCENARIO}.fastder.boundary_within_5bp"] == pytest.approx(75.0)
    assert found[f"sim.10M.{SCENARIO}.fastder.boundary_distance_median"] == pytest.approx(5)
    assert found[f"sim.10M.{SCENARIO}.fastder.locus_recall"] == pytest.approx(50.0)
    assert found[f"sim.10M.{SCENARIO}.fastder.stranded"] == pytest.approx(45.0)
    assert found[f"sim.10M.{SCENARIO}.fastder.strand_concordant"] == pytest.approx(75.0)


def test_runtime_takes_the_median_per_job_then_over_jobs(tree):
    found = values(tree)
    assert found["sim.10M.fastder.wall_median"] == pytest.approx(10.0)
    assert found["sim.10M.fastder.wall_max"] == pytest.approx(12.0)
    assert found["sim.10M.fastder.rss_median"] == pytest.approx(110.0)
    assert found["sim.10M.derfinder.wall_ratio_to_fastder"] == pytest.approx(10.0)
    assert found["sim.10M.scaling.cores4.wall"] == pytest.approx(6.0)
    assert found["sim.10M.scaling.cores4.rss"] == pytest.approx(2.0)


def test_an_unsampled_job_has_no_memory(tree):
    jobs = {job["job"]: job for job in benchmark_medians(str(tree[1] / "config_full_simulation"),
                                                         "run_fastder")}
    assert jobs["b"]["rss_mb"] is None


def test_replicates_get_their_own_names(tree):
    results, bench = tree
    make_run(results, "config_full_simulation_rep2")
    make_benchmark(bench, "config_full_simulation_rep2", "run_fastder", "a", [(9, 100)])
    numbers = collect(str(results), str(bench),
                      ["config_full_simulation", "config_full_simulation_rep2"])
    names = {n["name"] for n in numbers}
    assert f"sim.10M.rep2.{SCENARIO}.fastder.exon_prec" in names


def test_a_listed_config_with_a_missing_file_is_an_error(tree):
    results, bench = tree
    (results / "config_full_simulation" / "fuzzy_strand.csv").unlink()
    with pytest.raises(FileNotFoundError, match="fuzzy_strand.csv"):
        collect(str(results), str(bench), ["config_full_simulation"])
    with pytest.raises(FileNotFoundError):
        collect(str(results), str(bench), ["config_never_run"])


def test_threshold_ranges_and_comparison(tree):
    results, _ = tree
    write_rows(results / "ladder" / "threshold_range_summary.csv", [
        {"gene": "STMN2", "other_params": "coverage_tolerance=1.0", "thresholds_tested": 10,
         "thresholds_separating": 6, "lowest_separating": 0.005, "highest_separating": 0.2,
         "contiguous": 1},
        {"gene": "KCNQ2", "other_params": "coverage_tolerance=1.0", "thresholds_tested": 10,
         "thresholds_separating": 0, "lowest_separating": "", "highest_separating": "",
         "contiguous": ""},
    ])
    write_rows(results / "gtex" / "summary.csv", [
        {"tool": "fastder", "scenario": "BRAIN_1", "sample": "reference", "param_id": "mc1.0",
         "exon_prec": 50},
        {"tool": "fastder", "scenario": "HEART_1", "sample": "reference", "param_id": "mc1.0",
         "exon_prec": 54},
    ])
    found = values(tree, comparison_config="gtex", threshold_config="ladder")
    assert found["comparison.fastder.exon_prec"] == pytest.approx(52.0)
    assert found["threshold_range.STMN2.coverage_tolerance_1.0.highest_separating"] == 0.2
    assert found["threshold_range.KCNQ2.coverage_tolerance_1.0.thresholds_separating"] == 0
    assert "threshold_range.KCNQ2.coverage_tolerance_1.0.lowest_separating" not in found


def test_written_files(tree, tmp_path):
    results, bench = tree
    numbers = collect(str(results), str(bench), ["config_full_simulation"])
    write_csv(numbers, tmp_path / "n.csv")
    write_tex(numbers, tmp_path / "n.tex")
    with open(tmp_path / "n.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == len(numbers)
    assert set(rows[0]) == {"name", "value", "unit", "config", "source"}
    assert all(row["source"] and row["config"] for row in rows)
    tex = (tmp_path / "n.tex").read_text()
    assert tex.count(r"\@namedef{reported@") == len(numbers)
    assert rf"\@namedef{{reported@sim.10M.{SCENARIO}.fastder.exon_prec}}{{55.0}}" in tex


def test_tex_rounds_by_unit():
    assert tex_value({"value": 0.37149, "unit": "jaccard"}) == "0.371"
    assert tex_value({"value": 71.46, "unit": "percent"}) == "71.5"
    assert tex_value({"value": 4480.7, "unit": "MB"}) == "4481"


def test_ablation_junction_filter_and_annotation(tree):
    results, _ = tree
    sweep = results / "config_min_junction_reads_sweep"
    write_rows(sweep / "summary.csv", [
        {"tool": "fastder", "scenario": SCENARIO, "sample": "es",
         "param_id": f"mc0.005_ml10_pt5_mjr{v}_ns0", "exon_prec": 60 + v, "exon_sens": 40}
        for v in (0, 5)])
    write_rows(results / "config_unannotated_alignment" / "summary.csv", [
        {"tool": "fastder", "scenario": SCENARIO, "sample": "es", "param_id": FASTDER,
         "exon_prec": 48, "exon_sens": 35}])
    found = values(tree, junction_filter_config="config_min_junction_reads_sweep",
                   unannotated_config="config_unannotated_alignment")
    assert found[f"ablation.sim.10M.{SCENARIO}.stitched.exon_prec"] == pytest.approx(55.0)
    assert found[f"ablation.sim.10M.{SCENARIO}.unstitched.exon_prec"] == pytest.approx(1.0)
    assert found[f"junction_filter.{SCENARIO}.mjr5.exon_prec"] == pytest.approx(65.0)
    assert found[f"annotation.{SCENARIO}.unannotated.exon_prec"] == pytest.approx(48.0)
    assert found[f"annotation.{SCENARIO}.annotated.exon_prec"] == pytest.approx(55.0)


def test_a_listed_sweep_that_never_ran_is_an_error(tree):
    with pytest.raises(FileNotFoundError, match="config_unannotated_alignment"):
        values(tree, unannotated_config="config_unannotated_alignment")
