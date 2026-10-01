"""Coverage threshold with the best exon-level F1 on a threshold ladder.

Reads the summary.csv of a ladder run, in which each tool ran once per
min_coverage. F1 is the harmonic mean of gffcompare sensitivity and precision,
computed per tool, scenario and sample. The chosen threshold has the highest
exon-level F1 averaged over samples, scenarios and tools, so it favours no
single tool. A tie goes to the higher threshold.

Usage:
    python choose_threshold.py --summary <csv> --out-ladder <csv> --out-choice <csv>
"""
import argparse
import csv
from collections import defaultdict

from param_grid import parse_param_id

LEVELS = ("exon", "base")


def f1(sensitivity, precision):
    if sensitivity + precision == 0:
        return 0.0
    return 2 * sensitivity * precision / (sensitivity + precision)


def number(raw):
    """A gffcompare percentage, or None where it reports none."""
    if raw in (None, "", "nan", "-nan"):
        return None
    return float(raw)


def read_ladder(summary_csv):
    """One row per tool, scenario, sample and threshold, with F1 per level."""
    rows = []
    with open(summary_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            threshold = parse_param_id(row["param_id"]).get("min_coverage")
            if threshold is None:
                continue
            out = {"tool": row["tool"], "scenario": row["scenario"], "sample": row["sample"],
                   "min_coverage": float(threshold)}
            for level in LEVELS:
                sensitivity = number(row.get(f"{level}_sens"))
                precision = number(row.get(f"{level}_prec"))
                if sensitivity is None or precision is None:
                    break
                out[f"{level}_sens"] = sensitivity
                out[f"{level}_prec"] = precision
                out[f"{level}_f1"] = f1(sensitivity, precision)
            else:
                rows.append(out)
    return rows


def mean(values):
    values = list(values)
    return sum(values) / len(values)


def ladder_table(rows):
    """Mean over samples, per tool, scenario and threshold."""
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["tool"], row["scenario"], row["min_coverage"])].append(row)
    metrics = [f"{level}_{kind}" for level in LEVELS for kind in ("sens", "prec", "f1")]
    return [{"tool": tool, "scenario": scenario, "min_coverage": threshold, "samples": len(group),
             **{metric: mean(r[metric] for r in group) for metric in metrics}}
            for (tool, scenario, threshold), group in sorted(grouped.items())]


def choice_table(table):
    """Mean F1 per threshold over tools and scenarios, the best one flagged."""
    grouped = defaultdict(list)
    for row in table:
        grouped[row["min_coverage"]].append(row)
    choice = [{"min_coverage": threshold,
               "exon_f1": mean(r["exon_f1"] for r in group),
               "base_f1": mean(r["base_f1"] for r in group),
               "tools": len({r["tool"] for r in group})}
              for threshold, group in sorted(grouped.items())]
    if choice:
        best = max(choice, key=lambda r: (r["exon_f1"], r["min_coverage"]))
        for row in choice:
            row["chosen"] = int(row is best)
    return choice


def write_csv(rows, path, columns):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--out-ladder", required=True)
    parser.add_argument("--out-choice", required=True)
    args = parser.parse_args()

    table = ladder_table(read_ladder(args.summary))
    choice = choice_table(table)
    write_csv(table, args.out_ladder,
              ["tool", "scenario", "min_coverage", "samples",
               "exon_sens", "exon_prec", "exon_f1", "base_sens", "base_prec", "base_f1"])
    write_csv(choice, args.out_choice, ["min_coverage", "exon_f1", "base_f1", "tools", "chosen"])
    chosen = [row["min_coverage"] for row in choice if row["chosen"]]
    print(f"chosen min_coverage: {chosen[0] if chosen else 'none'}")


if __name__ == "__main__":
    main()
