"""Tidy table for the core scaling sweep.

run_fastder_scaling writes one benchmark TSV per core count. Memory is
reported next to wall time: each parsing thread holds one sample, so cores are
traded against memory.

Usage:
    python collect_scaling.py --bench-dir <dir> --out <csv>
"""
import argparse
import csv
import os.path as op
import re
from glob import glob

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
        for record in read_benchmark(path):
            wall = record.get("s")
            rss = record.get("max_rss")
            if wall in (None, ""):
                continue
            rows.append({
                "cores": cores,
                "wall_s": float(wall),
                # max_rss is in MB and is empty for a run that finished inside
                # the first sampling interval.
                "peak_rss_gb": float(rss) / MB_PER_GB if rss not in (None, "", "-") else "",
            })
    rows.sort(key=lambda r: r["cores"])
    return rows


def add_speedup(rows):
    """Speedup against the single-core point, when there is one."""
    single = [r["wall_s"] for r in rows if r["cores"] == 1]
    baseline = min(single) if single else None
    for row in rows:
        row["speedup"] = baseline / row["wall_s"] if baseline and row["wall_s"] else ""
    return rows


def write_csv(rows, path):
    columns = ["cores", "wall_s", "peak_rss_gb", "speedup"]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bench-dir", required=True,
                        help="benchmark directory of the config that ran the sweep")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rows = add_speedup(collect(args.bench_dir))
    write_csv(rows, args.out)
    print(f"wrote {len(rows)} rows to {args.out}")


if __name__ == "__main__":
    main()
