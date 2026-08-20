"""Tidy tables for the two single-axis fastder sweeps of the revision.

Both hold every fastder parameter at its shipped default and move one axis: the
ablation moves --no-stitch, the read-support sweep moves --min-junction-reads.
The panels plot these tables rather than recomputing, so a figure and the
numbers on disk cannot disagree.

Exon sensitivity and precision come from summary.csv, and the share of
predicted exon boundaries within 5 bp of a reference boundary from
fuzzy_distances.csv. Values are averaged over the samples of a scenario.

Usage:
    python collect_param_sweeps.py --axis no_stitch --results-root <dir> --out <csv>
    python collect_param_sweeps.py --axis min_junction_reads --results-root <dir> --out <csv>
"""
import argparse
import csv
import os
import os.path as op
import re
from collections import defaultdict

from param_grid import parse_param_id

# fastder's shipped defaults, from cpp/main.cpp. A sweep row is only comparable
# with the others if every axis but the swept one sits here. min_junction_reads
# is listed because a config may sweep both axes: without it, the ablation
# would average filtered and unfiltered runs into the same arm.
DEFAULTS = {"min_coverage": 0.05, "min_length": 10, "position_tolerance": 5,
            "min_junction_reads": 0}

BOUNDARY_WINDOW_BP = 5


def depth_of(run_dir):
    """Reads per sample in millions, from the results directory name.

    config_full_simulation is the 10M point and carries no suffix.
    """
    match = re.search(r"_([0-9]+)M$", op.basename(run_dir))
    return int(match.group(1)) if match else 10


def simulation_run_dirs(results_root, prefix="config_full_simulation"):
    if not op.isdir(results_root):
        return []
    return sorted(op.join(results_root, name) for name in os.listdir(results_root)
                  if name.startswith(prefix) and op.isdir(op.join(results_root, name)))


def comparable(combo, axis):
    """True when every parameter but the swept axis is at its default.

    --no-stitch makes position_tolerance inert, so the identifier of an
    unstitched run does not carry it: absence is accepted, a present but
    different value is not.
    """
    for name, default in DEFAULTS.items():
        if name == axis:
            continue
        if name in combo and combo[name] != default:
            return False
    return True


def axis_value(combo, axis):
    """The swept value, with the meaning an absent parameter has.

    An identifier written before --min-junction-reads existed describes a run
    with no junction filter, which is 0. An identifier with no ns component
    describes a stitched run.
    """
    if axis in combo:
        return combo[axis]
    return False if axis == "no_stitch" else 0


def read_rows(path):
    if not op.exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def accuracy_rows(run_dir, axis, tool):
    """Exon sensitivity and precision per swept value, from summary.csv."""
    out = defaultdict(list)
    for row in read_rows(op.join(run_dir, "summary.csv")):
        if row.get("tool") != tool:
            continue
        combo = parse_param_id(row.get("param_id", ""))
        if not comparable(combo, axis):
            continue
        value = axis_value(combo, axis)
        for metric, column in (("exon_sens", "exon_sens"), ("exon_prec", "exon_prec")):
            raw = row.get(column)
            if raw not in (None, ""):
                out[(row["scenario"], value, metric)].append(float(raw))
    return out


def boundary_rows(run_dir, axis, tool):
    """Share of boundaries within the window, from fuzzy_distances.csv."""
    hits = defaultdict(lambda: [0, 0])
    for row in read_rows(op.join(run_dir, "fuzzy_distances.csv")):
        if row.get("tool") != tool:
            continue
        combo = parse_param_id(row.get("param_id", ""))
        if not comparable(combo, axis):
            continue
        raw = row.get("distance")
        if raw in (None, ""):
            continue
        key = (row["scenario"], axis_value(combo, axis), "boundary_within_5bp")
        counts = hits[key]
        counts[1] += 1
        if abs(int(raw)) <= BOUNDARY_WINDOW_BP:
            counts[0] += 1
    return {key: [100.0 * hit / total] for key, (hit, total) in hits.items() if total}


def collect(results_root, axis, tool="fastder", prefix="config_full_simulation"):
    """One row per depth, scenario, swept value and metric."""
    rows = []
    for run_dir in simulation_run_dirs(results_root, prefix):
        depth = depth_of(run_dir)
        gathered = accuracy_rows(run_dir, axis, tool)
        for key, values in boundary_rows(run_dir, axis, tool).items():
            gathered[key].extend(values)
        ordered = sorted(gathered.items(), key=lambda kv: (kv[0][0], float(kv[0][1]), kv[0][2]))
        for (scenario, value, metric), values in ordered:
            rows.append({
                "depth_M": depth,
                "scenario": scenario,
                "tool": tool,
                axis: int(value) if isinstance(value, bool) else value,
                "metric": metric,
                "value": sum(values) / len(values),
                "n": len(values),
            })
    return rows


def write_csv(rows, path, axis):
    columns = ["depth_M", "scenario", "tool", axis, "metric", "value", "n"]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--axis", required=True,
                        choices=["no_stitch", "min_junction_reads"])
    parser.add_argument("--results-root", required=True,
                        help="workflow/results, holding one directory per config")
    parser.add_argument("--config-prefix", default="config_full_simulation",
                        help="results directories to read, by name prefix")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rows = collect(args.results_root, args.axis, prefix=args.config_prefix)
    write_csv(rows, args.out, args.axis)
    print(f"wrote {len(rows)} rows to {args.out}")


if __name__ == "__main__":
    main()
