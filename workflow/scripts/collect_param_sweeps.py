"""Tidy tables for the single-axis fastder comparisons.

Every parameter stays at its default while one axis moves: --no-stitch
for the ablation, --min-junction-reads for the read-support sweep, and the
STAR index, with or without the annotation, for annotated_index. Accuracy
comes from summary.csv, boundary distances from fuzzy_distances.csv, averaged
over the samples of a scenario.

The ablation has a third configuration, tool fastder_split: the stitched run
with every exon as its own record. It is graded per sample under
<run>/fastder_split and appears in neither summary file.

Usage:
    python collect_param_sweeps.py --axis no_stitch --results-root <dir> --out <csv>
    python collect_param_sweeps.py --axis min_junction_reads --results-root <dir> --out <csv>
    python collect_param_sweeps.py --axis annotated_index --results-root <dir> --out <csv>
"""
import argparse
import csv
import functools
import os
import os.path as op
import re
from collections import Counter, defaultdict
from glob import glob

from param_grid import comparable, parse_param_id
from parse_gffcompare import parse as parse_gffcompare_stats

BOUNDARY_WINDOW_BP = 5

SPLIT_TOOL = "fastder_split"

# The two alignments of the 10M reads, by whether the index held the annotation.
ANNOTATION_RUNS = ((1, "config_full_simulation"), (0, "config_unannotated_alignment"))
# The tools run on both alignments.
ANNOTATION_TOOLS = ("fastder", "derfinder", "megadepth_baseline")


def depth_of(run_dir):
    """Reads per sample in millions, from the results directory name.

    config_full_simulation is the 10M point and carries no suffix; neither do
    its replicates.
    """
    match = re.search(r"_([0-9]+)M$", op.basename(run_dir))
    return int(match.group(1)) if match else 10


def replicate_of(run_dir):
    """Draw of the simulation, from the _rep<N> suffix. An unsuffixed run is 1."""
    match = re.search(r"_rep([0-9]+)$", op.basename(run_dir))
    return int(match.group(1)) if match else 1


def simulation_run_dirs(results_root, prefix="config_full_simulation"):
    if not op.isdir(results_root):
        return []
    return sorted(op.join(results_root, name) for name in os.listdir(results_root)
                  if name.startswith(prefix) and op.isdir(op.join(results_root, name)))


def axis_value(combo, axis):
    """The swept value. An absent parameter means the published behaviour:
    no junction filter, stitching on."""
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


@functools.lru_cache(maxsize=None)
def boundary_distances(run_dir):
    """Absolute boundary distances of a run, counted per tool, scenario, sample
    and param_id. fuzzy_distances.csv has one row per boundary and runs to
    gigabytes, so it is read once per run and kept as counts."""
    path = op.join(run_dir, "fuzzy_distances.csv")
    counts = defaultdict(Counter)
    if not op.exists(path):
        return counts
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        if header is None:
            return counts
        tool, scenario, sample, param_id, distance = (
            header.index(name) for name in ("tool", "scenario", "sample", "param_id", "distance"))
        for row in reader:
            if row[distance] != "":
                counts[(row[tool], row[scenario], row[sample], row[param_id])][abs(int(row[distance]))] += 1
    return counts


def within_window_share(distance_counts):
    """Percent of the counted distances inside the window."""
    total = sum(distance_counts.values())
    hits = sum(n for distance, n in distance_counts.items() if distance <= BOUNDARY_WINDOW_BP)
    return 100.0 * hits / total


def boundary_rows(run_dir, axis, tool):
    """Share of boundaries within the window, from fuzzy_distances.csv."""
    pooled = defaultdict(Counter)
    for (row_tool, scenario, _, param_id), distance_counts in boundary_distances(run_dir).items():
        if row_tool != tool:
            continue
        combo = parse_param_id(param_id)
        if not comparable(combo, axis):
            continue
        pooled[(scenario, axis_value(combo, axis), "boundary_within_5bp")].update(distance_counts)
    return {key: [within_window_share(distance_counts)] for key, distance_counts in pooled.items()}


def within_window(distance):
    return abs(int(distance)) <= BOUNDARY_WINDOW_BP


def split_chain_rows(run_dir):
    """Exon accuracy and boundary share of the split-chain runs, per scenario."""
    out = defaultdict(list)
    hits = defaultdict(lambda: [0, 0])
    graded = op.join(run_dir, SPLIT_TOOL, "*", "*", "*")
    for stats in sorted(glob(op.join(graded, "gffcompare.stats"))):
        scenario = stats.split(os.sep)[-4]
        parsed = parse_gffcompare_stats(stats)
        for metric in ("exon_sens", "exon_prec"):
            if metric in parsed:
                out[(scenario, False, metric)].append(parsed[metric])
    for distances in sorted(glob(op.join(graded, "fuzzy_distances.csv"))):
        scenario = distances.split(os.sep)[-4]
        counts = hits[(scenario, False, "boundary_within_5bp")]
        for row in read_rows(distances):
            if row.get("distance") in (None, ""):
                continue
            counts[1] += 1
            counts[0] += within_window(row["distance"])
    for key, (hit, total) in hits.items():
        if total:
            out[key] = [100.0 * hit / total]
    return out


def collect(results_root, axis, tool="fastder", prefix="config_full_simulation"):
    """One row per depth, replicate, scenario, swept value and metric."""
    rows = []
    for run_dir in simulation_run_dirs(results_root, prefix):
        depth = depth_of(run_dir)
        replicate = replicate_of(run_dir)
        gathered = accuracy_rows(run_dir, axis, tool)
        for key, values in boundary_rows(run_dir, axis, tool).items():
            gathered[key].extend(values)
        by_tool = [(tool, gathered)]
        if axis == "no_stitch":
            by_tool.append((SPLIT_TOOL, split_chain_rows(run_dir)))
        for name, found in by_tool:
            ordered = sorted(found.items(), key=lambda kv: (kv[0][0], float(kv[0][1]), kv[0][2]))
            for (scenario, value, metric), values in ordered:
                rows.append({
                    "depth_M": depth,
                    "replicate": replicate,
                    "scenario": scenario,
                    "tool": name,
                    axis: int(value) if isinstance(value, bool) else value,
                    "metric": metric,
                    "value": sum(values) / len(values),
                    "n": len(values),
                })
    return rows


def collect_annotation(results_root, tools=ANNOTATION_TOOLS, runs=ANNOTATION_RUNS):
    """One row per tool, alignment, scenario and metric, at the reference point,
    fastder stitched."""
    rows = []
    for tool in tools:
        for annotated, name in runs:
            run_dir = op.join(results_root, name)
            gathered = accuracy_rows(run_dir, "no_stitch", tool)
            for key, values in boundary_rows(run_dir, "no_stitch", tool).items():
                gathered[key].extend(values)
            for (scenario, unstitched, metric), values in sorted(gathered.items()):
                if unstitched:
                    continue
                rows.append({
                    "depth_M": depth_of(run_dir),
                    "replicate": replicate_of(run_dir),
                    "scenario": scenario,
                    "tool": tool,
                    "annotated_index": annotated,
                    "metric": metric,
                    "value": sum(values) / len(values),
                    "n": len(values),
                })
    return rows


def write_csv(rows, path, axis):
    columns = ["depth_M", "replicate", "scenario", "tool", axis, "metric", "value", "n"]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--axis", required=True,
                        choices=["no_stitch", "min_junction_reads", "annotated_index"])
    parser.add_argument("--results-root", required=True,
                        help="workflow/results, holding one directory per config")
    parser.add_argument("--config-prefix", default="config_full_simulation",
                        help="results directories to read, by name prefix")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    if args.axis == "annotated_index":
        rows = collect_annotation(args.results_root)
    else:
        rows = collect(args.results_root, args.axis, prefix=args.config_prefix)
    write_csv(rows, args.out, args.axis)
    print(f"wrote {len(rows)} rows to {args.out}")


if __name__ == "__main__":
    main()
