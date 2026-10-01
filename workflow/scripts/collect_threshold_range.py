"""Thresholds at which fastder separates two groups at a set of loci.

A locus is separated at a parameter combination when the case group has a
called exon overlapping it and the control group has none. The summary gives,
per locus and per setting of the other parameters, the lowest and highest
min_coverage that separate the groups, and whether every threshold in between
does too.

Usage:
    python collect_threshold_range.py --loci <tsv> --case <dir> --control <dir> \
        --out <csv> --out-summary <csv>

Each group directory holds one <param_id>/output.gtf per parameter combination.
"""
import argparse
import csv
import os
import os.path as op
from collections import defaultdict

from param_grid import parse_param_id

SWEPT = "min_coverage"


def read_loci(path):
    """Loci with 1-based inclusive ce_start and ce_end, as GTF coordinates are."""
    with open(path, newline="") as fh:
        return [{"gene": r["gene"], "chrom": r["chrom"],
                 "start": int(r["ce_start"]), "end": int(r["ce_end"])}
                for r in csv.DictReader(fh, delimiter="\t")]


def read_exons(gtf_path):
    """Exon intervals per chromosome."""
    exons = defaultdict(list)
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) >= 5 and cols[2] == "exon":
                exons[cols[0]].append((int(cols[3]), int(cols[4])))
    return exons


def called(locus, exons):
    return any(start <= locus["end"] and end >= locus["start"]
               for start, end in exons.get(locus["chrom"], ()))


def param_ids(group_dir):
    return sorted(name for name in os.listdir(group_dir)
                  if op.exists(op.join(group_dir, name, "output.gtf")))


def other_params(combo):
    """The parameters held still, as text, so runs differing in them stay apart."""
    return ";".join(f"{name}={value}" for name, value in sorted(combo.items())
                    if name != SWEPT)


def collect(loci, case_dir, control_dir):
    """One row per locus and parameter combination present in both groups."""
    rows = []
    for param_id in param_ids(case_dir):
        control_gtf = op.join(control_dir, param_id, "output.gtf")
        if not op.exists(control_gtf):
            continue
        combo = parse_param_id(param_id)
        case_exons = read_exons(op.join(case_dir, param_id, "output.gtf"))
        control_exons = read_exons(control_gtf)
        for locus in loci:
            in_case = called(locus, case_exons)
            in_control = called(locus, control_exons)
            rows.append({
                "gene": locus["gene"],
                "param_id": param_id,
                SWEPT: combo.get(SWEPT, ""),
                "other_params": other_params(combo),
                "called_case": int(in_case),
                "called_control": int(in_control),
                "separates": int(in_case and not in_control),
            })
    return rows


def summarise(rows):
    """Lowest and highest separating threshold per locus and other parameters."""
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["gene"], row["other_params"])].append(row)
    summary = []
    for (gene, others), group in sorted(grouped.items()):
        group.sort(key=lambda r: float(r[SWEPT]))
        separating = [i for i, r in enumerate(group) if r["separates"]]
        if separating:
            lowest, highest = separating[0], separating[-1]
            contiguous = int(len(separating) == highest - lowest + 1)
            low_value, high_value = group[lowest][SWEPT], group[highest][SWEPT]
        else:
            contiguous, low_value, high_value = "", "", ""
        summary.append({
            "gene": gene,
            "other_params": others,
            "thresholds_tested": len(group),
            "thresholds_separating": len(separating),
            "lowest_separating": low_value,
            "highest_separating": high_value,
            "contiguous": contiguous,
        })
    return summary


def write_csv(rows, path, columns):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--loci", required=True)
    parser.add_argument("--case", required=True, help="tool output directory of the case group")
    parser.add_argument("--control", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--out-summary", required=True)
    args = parser.parse_args()

    rows = collect(read_loci(args.loci), args.case, args.control)
    write_csv(rows, args.out,
              ["gene", "param_id", SWEPT, "other_params",
               "called_case", "called_control", "separates"])
    write_csv(summarise(rows), args.out_summary,
              ["gene", "other_params", "thresholds_tested", "thresholds_separating",
               "lowest_separating", "highest_separating", "contiguous"])
    print(f"wrote {len(rows)} rows to {args.out}")


if __name__ == "__main__":
    main()
