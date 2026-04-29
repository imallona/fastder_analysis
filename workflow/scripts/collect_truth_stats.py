"""Aggregate per-transcript truth statistics from each scenario's splicing_variants.gff3 (the file used as the gffcompare reference). One row per truth transcript with scenario, sample, transcript_id, n_exons, total_exon_length, strand."""
import argparse
import csv
import os.path as op
import re
from collections import defaultdict


GFF_ATTR_RE = re.compile(r'(\w+)=([^;]+)')


def parse_gff(gff_path, scenario, sample):
    transcripts_n_exons = defaultdict(int)
    transcripts_length = defaultdict(int)
    transcripts_strand = {}
    transcripts_chrom = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            attrs = dict(GFF_ATTR_RE.findall(cols[8]))
            tx = attrs.get("transcript_id", "")
            if not tx:
                continue
            transcripts_n_exons[tx] += 1
            transcripts_length[tx] += int(cols[4]) - int(cols[3]) + 1
            transcripts_strand[tx] = cols[6]
            transcripts_chrom[tx] = cols[0]
    rows = []
    for tx, n in transcripts_n_exons.items():
        rows.append({
            "scenario": scenario,
            "sample": sample,
            "transcript_id": tx,
            "chrom": transcripts_chrom[tx],
            "n_exons": n,
            "total_exon_length": transcripts_length[tx],
            "strand": transcripts_strand[tx],
        })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", action="append", required=True,
                    help="Path to a splicing_variants.gff3 (repeat per scenario+sample)")
    ap.add_argument("--scenario", action="append", required=True,
                    help="scenario tag for each --gff (same order)")
    ap.add_argument("--sample", action="append", required=True,
                    help="sample tag for each --gff (same order)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    if not (len(args.gff) == len(args.scenario) == len(args.sample)):
        ap.error("--gff, --scenario, and --sample must be paired")

    fieldnames = ["scenario", "sample", "transcript_id", "chrom", "n_exons",
                  "total_exon_length", "strand"]
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for gff, scenario, sample in zip(args.gff, args.scenario, args.sample):
            if not op.isfile(gff):
                continue
            for row in parse_gff(gff, scenario, sample):
                w.writerow(row)


if __name__ == "__main__":
    main()
