"""Every value the text quotes, in one table.

Reads the result and benchmark files of the listed configs and writes one row
per value: name, value, unit, config and source file. A TeX file defines one
macro per row, read as \\reported{<name>}. A config that is listed and has a
file missing is an error; nothing is written blank.

Tools are read at the reference point of param_grid.py: fastder stitched,
derfinder and the megadepth baseline at the same coverage threshold,
groHMM at the combination with the highest median Jaccard in that config.
Values are means over the samples of a scenario unless the name says median.

Usage:
    python collect_reported_numbers.py --results-root <dir> --bench-root <dir> \\
        --simulation <config> [<config> ...] [--comparison <config>] \\
        [--runtime <config> ...] [--threshold-range <config>] \\
        [--junction-filter <config>] [--unannotated <config>] \\
        --out-csv <csv> --out-tex <tex>
"""
import argparse
import csv
import os.path as op
import re
from collections import defaultdict
from glob import glob
from statistics import mean, median

import collect_param_sweeps as sweeps
from collect_param_sweeps import depth_of, replicate_of
from param_grid import REFERENCE, parse_param_id, stitched_at_reference

TOOLS = ["fastder", "derfinder", "megadepth_baseline", "grohmm"]
TIMED_RULES = {tool: f"run_{tool}" for tool in TOOLS}
LEVELS = ["base", "exon", "intron", "intron_chain", "transcript", "locus"]
BOUNDARY_WINDOW_BP = 5
LOCUS_RECALL_THRESHOLD = 0.5
MB_PER_GB = 1024.0

# Decimals in the TeX file, by unit. The CSV keeps full precision.
DECIMALS = {"percent": 1, "s": 1, "MB": 0, "GB": 1, "jaccard": 3, "bp": 0,
            "ratio": 1, "CPM": 3, "count": 0}


def read_rows(path, delimiter=","):
    if not op.exists(path):
        raise FileNotFoundError(f"expected input is missing: {path}")
    with open(path, newline="") as fh:
        yield from csv.DictReader(fh, delimiter=delimiter)


def number(name, value, unit, config, source):
    return {"name": name, "value": value, "unit": unit, "config": config,
            "source": op.basename(source)}


def run_label(config):
    """sim.<depth>M, with .rep<N> for a replicate after the first."""
    replicate = replicate_of(config)
    suffix = f".rep{replicate}" if replicate > 1 else ""
    return f"sim.{depth_of(config)}M{suffix}"


def reference_params(run_dir):
    """Parameter identifier per tool, as the text compares them."""
    fastder, grohmm_jaccard = None, defaultdict(list)
    for row in read_rows(op.join(run_dir, "fuzzy_jaccard.csv")):
        tool, param_id = row["tool"], row["param_id"]
        if tool == "fastder" and fastder is None and stitched_at_reference(parse_param_id(param_id)):
            fastder = param_id
        elif tool == "grohmm":
            grohmm_jaccard[param_id].append(float(row["jaccard"]))
    coverage = REFERENCE["min_coverage"]
    params = {
        "fastder": fastder,
        "derfinder": f"mc{coverage}_pt{REFERENCE['position_tolerance']}",
        "megadepth_baseline": f"mc{coverage}",
    }
    if grohmm_jaccard:
        params["grohmm"] = max(sorted(grohmm_jaccard), key=lambda p: median(grohmm_jaccard[p]))
    return {tool: param_id for tool, param_id in params.items() if param_id}


def at_reference(rows, params):
    for row in rows:
        if params.get(row["tool"]) == row["param_id"]:
            yield row


def accuracy(config, run_dir, params):
    """gffcompare sensitivity and precision, per scenario and tool, and per
    sample for fastder."""
    source = op.join(run_dir, "summary.csv")
    values = defaultdict(list)
    per_sample = {}
    for row in at_reference(read_rows(source), params):
        for level in LEVELS:
            for side in ("sens", "prec"):
                raw = row.get(f"{level}_{side}")
                if raw in (None, "", "nan", "-nan"):
                    continue
                key = (row["scenario"], row["tool"], f"{level}_{side}")
                values[key].append(float(raw))
                if row["tool"] == "fastder" and level == "exon":
                    per_sample[(row["scenario"], row["sample"], f"{level}_{side}")] = float(raw)
    label = run_label(config)
    out = [number(f"{label}.{scenario}.{tool}.{metric}", mean(found), "percent", config, source)
           for (scenario, tool, metric), found in sorted(values.items())]
    out += [number(f"{label}.{scenario}.fastder.{sample}.{metric}", value, "percent", config, source)
            for (scenario, sample, metric), value in sorted(per_sample.items())]
    return out


def overlap(config, run_dir, params):
    """Median exonic Jaccard of a reference transcript's best call."""
    source = op.join(run_dir, "fuzzy_jaccard.csv")
    values = defaultdict(list)
    for row in at_reference(read_rows(source), params):
        values[(row["scenario"], row["tool"])].append(float(row["jaccard"]))
    label = run_label(config)
    return [number(f"{label}.{scenario}.{tool}.jaccard_median", median(found), "jaccard",
                   config, source)
            for (scenario, tool), found in sorted(values.items())]


def boundaries(config, run_dir, params):
    """Share of boundaries within the window, averaged over samples, and the
    median absolute distance."""
    source = op.join(run_dir, "fuzzy_distances.csv")
    hits = defaultdict(lambda: [0, 0])
    distances = defaultdict(list)
    for row in at_reference(read_rows(source), params):
        if row["distance"] in (None, ""):
            continue
        distance = abs(int(row["distance"]))
        counts = hits[(row["scenario"], row["tool"], row["sample"])]
        counts[1] += 1
        counts[0] += distance <= BOUNDARY_WINDOW_BP
        distances[(row["scenario"], row["tool"])].append(distance)
    shares = defaultdict(list)
    for (scenario, tool, _), (hit, total) in hits.items():
        shares[(scenario, tool)].append(100.0 * hit / total)
    label = run_label(config)
    out = []
    for (scenario, tool), found in sorted(shares.items()):
        out.append(number(f"{label}.{scenario}.{tool}.boundary_within_{BOUNDARY_WINDOW_BP}bp",
                          mean(found), "percent", config, source))
        out.append(number(f"{label}.{scenario}.{tool}.boundary_distance_median",
                          median(distances[(scenario, tool)]), "bp", config, source))
    return out


def locus_recall(config, run_dir, params):
    source = op.join(run_dir, "fuzzy_locus_recall.csv")
    values = defaultdict(list)
    for row in at_reference(read_rows(source), params):
        if abs(float(row["threshold"]) - LOCUS_RECALL_THRESHOLD) < 1e-9:
            values[(row["scenario"], row["tool"])].append(100.0 * float(row["recall"]))
    label = run_label(config)
    return [number(f"{label}.{scenario}.{tool}.locus_recall", mean(found), "percent",
                   config, source)
            for (scenario, tool), found in sorted(values.items())]


def strand(config, run_dir, params):
    """Share of fastder calls that carry a strand, and of those matched to a
    reference transcript, the share on its strand."""
    source = op.join(run_dir, "fuzzy_strand.csv")
    counts = defaultdict(lambda: defaultdict(int))
    for row in at_reference(read_rows(source), params):
        if row["tool"] == "fastder":
            counts[row["scenario"]][row["category"]] += int(row["n_fastder_transcripts"])
    label = run_label(config)
    out = []
    for scenario, found in sorted(counts.items()):
        total = sum(found.values())
        matched = found["concordant"] + found["discordant"]
        if total:
            out.append(number(f"{label}.{scenario}.fastder.stranded",
                              100.0 * (total - found["unstranded"]) / total,
                              "percent", config, source))
        if matched:
            out.append(number(f"{label}.{scenario}.fastder.strand_concordant",
                              100.0 * found["concordant"] / matched, "percent", config, source))
    return out


def benchmark_medians(bench_dir, rule):
    """Per job of a rule, the median wall time and peak memory over its repeats."""
    jobs = []
    for path in sorted(glob(op.join(bench_dir, rule, "*.tsv"))):
        rows = list(read_rows(path, delimiter="\t"))
        walls = [float(r["s"]) for r in rows if r.get("s") not in (None, "")]
        rss = [float(r["max_rss"]) for r in rows if r.get("max_rss") not in (None, "", "-")]
        if walls:
            jobs.append({"job": op.basename(path)[: -len(".tsv")], "wall_s": median(walls),
                         "rss_mb": median(rss) if rss else None})
    return jobs


def runtime(config, bench_dir, label):
    """Median, lowest and highest wall time per tool over its jobs, peak memory,
    and how many times longer each tool takes than fastder."""
    if not op.isdir(bench_dir):
        raise FileNotFoundError(f"expected input is missing: {bench_dir}")
    out, medians = [], {}
    for tool, rule in TIMED_RULES.items():
        jobs = benchmark_medians(bench_dir, rule)
        if not jobs:
            continue
        source = op.join(bench_dir, rule)
        walls = [job["wall_s"] for job in jobs]
        medians[tool] = median(walls)
        for stat, value in (("median", medians[tool]), ("min", min(walls)), ("max", max(walls))):
            out.append(number(f"{label}.{tool}.wall_{stat}", value, "s", config, source))
        rss = [job["rss_mb"] for job in jobs if job["rss_mb"] is not None]
        if rss:
            out.append(number(f"{label}.{tool}.rss_median", median(rss), "MB", config, source))
            out.append(number(f"{label}.{tool}.rss_max", max(rss), "MB", config, source))
    for tool, value in medians.items():
        if tool != "fastder" and medians.get("fastder"):
            out.append(number(f"{label}.{tool}.wall_ratio_to_fastder", value / medians["fastder"],
                              "ratio", config, op.join(bench_dir, TIMED_RULES[tool])))
    for job in benchmark_medians(bench_dir, "run_fastder_scaling"):
        cores = re.sub(r"\D", "", job["job"])
        source = op.join(bench_dir, "run_fastder_scaling")
        out.append(number(f"{label}.scaling.cores{cores}.wall", job["wall_s"], "s", config, source))
        if job["rss_mb"] is not None:
            out.append(number(f"{label}.scaling.cores{cores}.rss", job["rss_mb"] / MB_PER_GB,
                              "GB", config, source))
    return out


def comparison(config, run_dir):
    """gffcompare against the annotation, mean over the scenarios of a
    real-data comparison, each tool at its only parameter combination."""
    source = op.join(run_dir, "summary.csv")
    values = defaultdict(list)
    for row in read_rows(source):
        for level in LEVELS:
            for side in ("sens", "prec"):
                raw = row.get(f"{level}_{side}")
                if raw not in (None, "", "nan", "-nan"):
                    values[(row["tool"], f"{level}_{side}")].append(float(raw))
    return [number(f"comparison.{tool}.{metric}", mean(found), "percent", config, source)
            for (tool, metric), found in sorted(values.items())]


def threshold_ranges(config, run_dir):
    source = op.join(run_dir, "threshold_range_summary.csv")
    out = []
    for row in read_rows(source):
        setting = re.sub(r"[^A-Za-z0-9.]+", "_", row["other_params"]) or "defaults"
        name = f"threshold_range.{row['gene']}.{setting}"
        for column in ("lowest_separating", "highest_separating"):
            if row[column] != "":
                out.append(number(f"{name}.{column}", float(row[column]), "CPM", config, source))
        out.append(number(f"{name}.thresholds_separating", int(row["thresholds_separating"]),
                          "count", config, source))
    return out


ABLATION_NAMES = {("fastder", 0): "stitched", ("fastder", 1): "unstitched",
                  (sweeps.SPLIT_TOOL, 0): "split"}


def require_run(results_root, config):
    run_dir = op.join(results_root, config)
    if not op.isdir(run_dir):
        raise FileNotFoundError(f"expected input is missing: {run_dir}")
    return run_dir


def ablation(results_root, simulations):
    """Exon accuracy and boundary share of the three ablation configurations."""
    by_run = {(depth_of(config), replicate_of(config)): config for config in simulations}
    out = []
    for row in sweeps.collect(results_root, "no_stitch"):
        config = by_run.get((row["depth_M"], row["replicate"]))
        if config is None:
            continue
        arm = ABLATION_NAMES[(row["tool"], row["no_stitch"])]
        out.append(number(f"ablation.{run_label(config)}.{row['scenario']}.{arm}.{row['metric']}",
                          row["value"], "percent", config, "summary.csv"))
    return out


def junction_filter(results_root, config):
    require_run(results_root, config)
    return [number(f"junction_filter.{row['scenario']}.mjr{row['min_junction_reads']}.{row['metric']}",
                   row["value"], "percent", config, "summary.csv")
            for row in sweeps.collect(results_root, "min_junction_reads", prefix=config)]


def annotation(results_root, annotated_config, unannotated_config):
    """Each tool at the reference point on the annotated and the unannotated alignment."""
    require_run(results_root, unannotated_config)
    runs = ((1, annotated_config), (0, unannotated_config))
    return [number(f"annotation.{row['scenario']}.{row['tool']}."
                   f"{'annotated' if row['annotated_index'] else 'unannotated'}.{row['metric']}",
                   row["value"], "percent",
                   annotated_config if row["annotated_index"] else unannotated_config,
                   "summary.csv")
            for row in sweeps.collect_annotation(results_root, runs=runs)]


def collect(results_root, bench_root, simulations, comparison_config=None,
            runtime_configs=(), threshold_config=None, junction_filter_config=None,
            unannotated_config=None):
    numbers = []
    for config in simulations:
        run_dir = op.join(results_root, config)
        params = reference_params(run_dir)
        for section in (accuracy, overlap, boundaries, locus_recall, strand):
            numbers += section(config, run_dir, params)
        # Later replicates are run for accuracy only and are not timed.
        if replicate_of(config) == 1:
            numbers += runtime(config, op.join(bench_root, config), run_label(config))
    if comparison_config:
        numbers += comparison(comparison_config, op.join(results_root, comparison_config))
    for config in runtime_configs:
        numbers += runtime(config, op.join(bench_root, config), f"runtime.{config}")
    if threshold_config:
        numbers += threshold_ranges(threshold_config, op.join(results_root, threshold_config))
    numbers += ablation(results_root, simulations)
    if junction_filter_config:
        numbers += junction_filter(results_root, junction_filter_config)
    if unannotated_config:
        numbers += annotation(results_root, simulations[0], unannotated_config)
    names = [n["name"] for n in numbers]
    repeated = sorted({name for name in names if names.count(name) > 1})
    if repeated:
        raise ValueError(f"names are not unique: {repeated[:5]}")
    return numbers


def write_csv(numbers, path):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["name", "value", "unit", "config", "source"])
        writer.writeheader()
        writer.writerows(numbers)


def tex_value(row):
    if isinstance(row["value"], str):
        return row["value"]
    return f"{row['value']:.{DECIMALS[row['unit']]}f}"


def write_tex(numbers, path):
    """One macro per value. An unknown name stops the build."""
    lines = ["% Generated by workflow/scripts/collect_reported_numbers.py. Do not edit.",
             r"\makeatletter",
             r"\newcommand{\reported}[1]{\@ifundefined{reported@#1}"
             r"{\PackageError{reported}{unknown value #1}{}}{\@nameuse{reported@#1}}}"]
    lines += [f"\\@namedef{{reported@{row['name']}}}{{{tex_value(row)}}}" for row in numbers]
    lines.append(r"\makeatother")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results-root", required=True)
    parser.add_argument("--bench-root", required=True)
    parser.add_argument("--simulation", nargs="+", required=True)
    parser.add_argument("--comparison")
    parser.add_argument("--runtime", nargs="*", default=[])
    parser.add_argument("--threshold-range")
    parser.add_argument("--junction-filter")
    parser.add_argument("--unannotated")
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--out-tex", required=True)
    args = parser.parse_args()

    numbers = collect(args.results_root, args.bench_root, args.simulation,
                      args.comparison, args.runtime, args.threshold_range,
                      args.junction_filter, args.unannotated)
    write_csv(numbers, args.out_csv)
    write_tex(numbers, args.out_tex)
    print(f"wrote {len(numbers)} values to {args.out_csv} and {args.out_tex}")


if __name__ == "__main__":
    main()
