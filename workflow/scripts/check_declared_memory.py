#!/usr/bin/env python3
"""List rules whose measured peak memory exceeds their declared mem_mb.

Reads the benchmark files of one config and the mem_mb each rule declares.
Exits 1 if any rule used more than it declared.
"""

import argparse
import csv
import re
import sys
from pathlib import Path

RULE_HEADER = re.compile(r"^(?:rule|checkpoint)\s+(\w+):")
MEM_MB = re.compile(r"^\s*mem_mb=(\d+),")


def declared_memory(rules_dir):
    """Rule name to its declared mem_mb, for rules that declare a number."""
    declared = {}
    for path in sorted(Path(rules_dir).glob("*.smk")):
        rule = None
        for line in path.read_text().splitlines():
            header = RULE_HEADER.match(line)
            if header:
                rule = header.group(1)
                continue
            mem = MEM_MB.match(line)
            if rule and mem:
                declared[rule] = int(mem.group(1))
    return declared


def rule_of(benchmark, benchmark_dir, rules):
    """The rule a benchmark file belongs to: the longest rule name its path starts with."""
    relative = Path(benchmark).relative_to(benchmark_dir).as_posix()
    matches = [rule for rule in rules if relative.startswith(rule)]
    return max(matches, key=len) if matches else None


def measured_memory(benchmark_dir, rules):
    """Rule name to its highest max_rss in MB and the file that holds it."""
    measured = {}
    for benchmark in sorted(Path(benchmark_dir).rglob("*.tsv")):
        rule = rule_of(benchmark, benchmark_dir, rules)
        if rule is None:
            continue
        with open(benchmark, newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                try:
                    max_rss = float(row["max_rss"])
                except (KeyError, ValueError):
                    continue
                if max_rss > measured.get(rule, (0.0, None))[0]:
                    measured[rule] = (max_rss, str(benchmark))
    return measured


def excesses(declared, measured):
    """Rows of (rule, declared mem_mb, measured MB, file), largest excess first."""
    rows = [
        (rule, declared[rule], max_rss, benchmark)
        for rule, (max_rss, benchmark) in measured.items()
        if max_rss > declared[rule]
    ]
    return sorted(rows, key=lambda row: row[2] / row[1], reverse=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--rules", required=True, help="directory of rule files")
    parser.add_argument("--benchmarks", required=True, help="benchmark directory of one config")
    args = parser.parse_args()

    declared = declared_memory(args.rules)
    measured = measured_memory(args.benchmarks, declared)
    rows = excesses(declared, measured)
    print("rule\tdeclared_mb\tmeasured_mb\tbenchmark")
    for rule, declared_mb, measured_mb, benchmark in rows:
        print(f"{rule}\t{declared_mb}\t{measured_mb:.0f}\t{benchmark}")
    return 1 if rows else 0


if __name__ == "__main__":
    sys.exit(main())
