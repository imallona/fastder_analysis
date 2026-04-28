"""Aggregate per-transcript stats from each fastder GTF output into one CSV.

Snakemake calls this with one argument per parameter combination: a path to a
.gtf_path file produced by run_fastder. Each .gtf_path file contains a single
line, the absolute path of the FASTDER_RESULT_*.gtf for that param_id.
The output CSV has one row per transcript, with columns:
  param_id, chrom, transcript_id, strand, n_exons, total_exon_length, score
"""
import argparse
import csv
import os.path as op
import re

TRANSCRIPT_ID_RE = re.compile(r'transcript_id "([^"]+)"')


def parse_gtf(gtf_path, param_id):
    """Yield one dict per transcript in the gtf."""
    current = None
    with open(gtf_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            chrom, _, _, start, end, score, strand, _, attrs = parts
            m = TRANSCRIPT_ID_RE.search(attrs)
            tx_id = m.group(1) if m else ""
            if current is None or current["transcript_id"] != tx_id:
                if current is not None:
                    yield current
                current = {
                    "param_id": param_id,
                    "chrom": chrom,
                    "transcript_id": tx_id,
                    "strand": strand,
                    "n_exons": 0,
                    "total_exon_length": 0,
                    "score": float(score),
                }
            current["n_exons"] += 1
            current["total_exon_length"] += int(end) - int(start) + 1
        if current is not None:
            yield current


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gtf-path-file", action="append", required=True,
                   help="Path to a .gtf_path file (repeat for each param combo)")
    p.add_argument("--param-id", action="append", required=True,
                   help="param_id for each --gtf-path-file (in same order)")
    p.add_argument("--out", required=True, help="Output CSV path")
    args = p.parse_args()
    if len(args.gtf_path_file) != len(args.param_id):
        p.error("--gtf-path-file and --param-id must be paired")

    rows = []
    for gtf_path_file, param_id in zip(args.gtf_path_file, args.param_id):
        with open(gtf_path_file) as f:
            gtf_path = f.read().strip()
        if not op.isfile(gtf_path):
            raise FileNotFoundError(f"missing GTF for param {param_id}: {gtf_path}")
        rows.extend(parse_gtf(gtf_path, param_id))

    fieldnames = ["param_id", "chrom", "transcript_id", "strand",
                  "n_exons", "total_exon_length", "score"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
