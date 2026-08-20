"""The scaling table records how fastder uses the cores it is given.

The trap is the empty max_rss cell. Snakemake samples every 0.5 s for the first
thirty samples, so a run finishing before the first sample has a wall clock and
no memory figure. Reading that empty string as zero draws a memory curve
diving to the floor exactly where the run is fastest.
"""

import csv

import pytest

from collect_scaling import add_speedup, collect, cores_of, write_csv

BENCH_COLUMNS = ["s", "h:m:s", "max_rss", "max_vms", "max_uss", "max_pss",
                 "io_in", "io_out", "mean_load", "cpu_time"]


def write_benchmark(path, wall_s, max_rss):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=BENCH_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerow({"s": wall_s, "h:m:s": "0:00:10", "max_rss": max_rss,
                         "max_vms": "", "max_uss": "", "max_pss": "",
                         "io_in": "", "io_out": "", "mean_load": "", "cpu_time": ""})


def make_sweep(tmp_path, points):
    for cores, wall_s, max_rss in points:
        write_benchmark(tmp_path / "run_fastder_scaling" / f"cores{cores}.tsv",
                        wall_s, max_rss)
    return tmp_path


def test_cores_come_from_the_file_name():
    assert cores_of("/x/run_fastder_scaling/cores16.tsv") == 16
    assert cores_of("/x/run_fastder_scaling/other.tsv") is None


def test_rows_are_sorted_by_core_count(tmp_path):
    make_sweep(tmp_path, [(8, 20.0, 8192), (1, 100.0, 4096), (4, 30.0, 6144)])
    rows = collect(str(tmp_path))
    assert [r["cores"] for r in rows] == [1, 4, 8]


def test_memory_is_reported_in_gibibytes(tmp_path):
    make_sweep(tmp_path, [(1, 100.0, 15872)])
    row = collect(str(tmp_path))[0]
    assert row["peak_rss_gb"] == pytest.approx(15.5)


def test_an_unsampled_run_has_no_memory_rather_than_zero(tmp_path):
    make_sweep(tmp_path, [(16, 2.0, "")])
    row = collect(str(tmp_path))[0]
    assert row["peak_rss_gb"] == ""
    assert row["wall_s"] == pytest.approx(2.0)


def test_speedup_is_against_the_single_core_point(tmp_path):
    make_sweep(tmp_path, [(1, 100.0, 4096), (4, 25.0, 6144)])
    rows = add_speedup(collect(str(tmp_path)))
    by_cores = {r["cores"]: r for r in rows}
    assert by_cores[1]["speedup"] == pytest.approx(1.0)
    assert by_cores[4]["speedup"] == pytest.approx(4.0)


def test_speedup_is_blank_without_a_single_core_point(tmp_path):
    make_sweep(tmp_path, [(4, 25.0, 6144)])
    rows = add_speedup(collect(str(tmp_path)))
    assert rows[0]["speedup"] == ""


def test_written_csv_carries_every_column(tmp_path):
    make_sweep(tmp_path, [(1, 100.0, 4096)])
    out = tmp_path / "scaling.csv"
    write_csv(add_speedup(collect(str(tmp_path))), out)
    with open(out) as fh:
        header = next(csv.reader(fh))
    assert header == ["cores", "wall_s", "peak_rss_gb", "speedup"]


def test_missing_sweep_is_not_an_error(tmp_path):
    assert collect(str(tmp_path)) == []
