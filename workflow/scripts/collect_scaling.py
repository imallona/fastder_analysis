"""Tidy table for the core scaling sweep.

run_fastder_scaling writes one benchmark TSV per core count, with one row per
repeat; each core count is reported by its median. Memory is reported next to
wall time: each parsing thread holds one sample, so cores are traded against
memory.

A sweep saturates where its workload has nothing left to spread: parsing
at the number of samples, averaging at the number of chromosomes. Each
workload is given with both, so the figure can mark them.

Usage:
    python collect_scaling.py --out <csv> \
        --workload <label> <bench dir> <samples> <chromosomes> [--workload ...]
"""
import argparse
import csv
import os.path as op
import re
from glob import glob
from statistics import median

MB_PER_GB = 1024.0


def cores_of(path):
    match = re.search(r"cores([0-9]+)\.tsv$", op.basename(path))
    if not match:
        return None
    return int(match.group(1))


def read_benchmark(path):
    """Snakemake's benchmark TSV: one header row, one row per repeat."""
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows


def collect(bench_dir):
    rows = []
    for path in sorted(glob(op.join(bench_dir, "run_fastder_scaling", "cores*.tsv"))):
        cores = cores_of(path)
        if cores is None:
            continue
        records = read_benchmark(path)
        walls = [float(r["s"]) for r in records if r.get("s") not in (None, "")]
        # max_rss is in MB and is empty for a run that finished inside the
        # first sampling interval.
        rss = [float(r["max_rss"]) for r in records
               if r.get("max_rss") not in (None, "", "-")]
        if not walls:
            continue
        rows.append({
            "cores": cores,
            "wall_s": median(walls),
            "peak_rss_gb": median(rss) / MB_PER_GB if rss else "",
            "repeats": len(walls),
        })
    rows.sort(key=lambda r: r["cores"])
    return rows


def add_speedup(rows):
    """Speedup against the single-core point, when there is one."""
    baseline = next((r["wall_s"] for r in rows if r["cores"] == 1), None)
    for row in rows:
        row["speedup"] = baseline / row["wall_s"] if baseline and row["wall_s"] else ""
    return rows


def collect_workload(label, bench_dir, samples, chromosomes):
    rows = add_speedup(collect(bench_dir))
    for row in rows:
        row.update(workload=label, samples=int(samples), chromosomes=int(chromosomes))
    return rows


def write_csv(rows, path):
    columns = ["workload", "cores", "wall_s", "peak_rss_gb", "repeats", "speedup",
               "samples", "chromosomes"]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workload", nargs=4, action="append", required=True,
                        metavar=("LABEL", "BENCH_DIR", "SAMPLES", "CHROMOSOMES"),
                        help="one sweep: its label, the benchmark directory of the "
                             "config that ran it, and its sample and chromosome counts")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rows = [row for workload in args.workload for row in collect_workload(*workload)]
    write_csv(rows, args.out)
    print(f"wrote {len(rows)} rows to {args.out}")


if __name__ == "__main__":
    main()
