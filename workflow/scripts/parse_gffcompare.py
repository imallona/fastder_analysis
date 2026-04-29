"""Parse a gffcompare.stats file into a flat dict of named metrics.

The output dict is intended to be merged into one row of the summary CSV.
Keys cover all six "Level" rows (sensitivity and precision), the three
matching counts, the missed and novel counts (count, total, and percent),
and the header summary numbers.
"""
import argparse
import csv
import json
import re
import sys

LEVEL_RE = re.compile(r"^\s*([A-Za-z][A-Za-z ]+) level\s*:\s*([\d.]+)\s*\|\s*([\d.]+)")
MATCHING_RE = re.compile(r"^\s*Matching ([a-z ]+):\s*(\d+)\s*$")
MISSED_NOVEL_RE = re.compile(r"^\s*(Missed|Novel) ([a-z]+):\s*(\d+)/(\d+)\s*\(\s*([\d.]+)%\s*\)")
QUERY_HEADER_RE = re.compile(r"^#\s*Query mRNAs\s*:\s*(\d+)\s+in\s+(\d+)\s+loci\s*\(\s*(\d+)\s+multi-exon\s+transcripts\)")
REF_HEADER_RE = re.compile(r"^#\s*Reference mRNAs\s*:\s*(\d+)\s+in\s+(\d+)\s+loci\s*\(\s*(\d+)\s+multi-exon\)")
SUPER_LOCI_RE = re.compile(r"^#\s*Super-loci w/ reference transcripts:\s*(\d+)")


def slug(name):
    """Normalise a level name into a snake_case key fragment."""
    return name.strip().lower().replace(" ", "_")


def parse(path):
    out = {}
    with open(path) as f:
        for line in f:
            m = LEVEL_RE.match(line)
            if m:
                key = slug(m.group(1))
                out[f"{key}_sens"] = float(m.group(2))
                out[f"{key}_prec"] = float(m.group(3))
                continue
            m = MATCHING_RE.match(line)
            if m:
                key = slug(m.group(1))
                out[f"matching_{key}"] = int(m.group(2))
                continue
            m = MISSED_NOVEL_RE.match(line)
            if m:
                kind = m.group(1).lower()
                what = m.group(2).lower()
                out[f"{kind}_{what}_count"] = int(m.group(3))
                out[f"{kind}_{what}_total"] = int(m.group(4))
                out[f"{kind}_{what}_pct"]   = float(m.group(5))
                continue
            m = QUERY_HEADER_RE.match(line)
            if m:
                out["query_mrnas"] = int(m.group(1))
                out["query_loci"] = int(m.group(2))
                out["query_multi_exon"] = int(m.group(3))
                continue
            m = REF_HEADER_RE.match(line)
            if m:
                out["ref_mrnas"] = int(m.group(1))
                out["ref_loci"] = int(m.group(2))
                out["ref_multi_exon"] = int(m.group(3))
                continue
            m = SUPER_LOCI_RE.match(line)
            if m:
                out["super_loci"] = int(m.group(1))
                continue
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("path")
    p.add_argument("--json", action="store_true",
                   help="Emit JSON to stdout (default: one key=value per line)")
    args = p.parse_args()
    fields = parse(args.path)
    if args.json:
        json.dump(fields, sys.stdout, indent=2, sort_keys=True)
        sys.stdout.write("\n")
    else:
        for k, v in sorted(fields.items()):
            print(f"{k}\t{v}")


if __name__ == "__main__":
    main()
