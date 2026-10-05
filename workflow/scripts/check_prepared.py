#!/usr/bin/env python3
"""Check that a tree prepared elsewhere is ready for its timed pass.

--commit-file: the commit the tree was prepared at must be the one checked out.
--manifest: every file listed there must be under --root with the same size.
--plan: a Snakemake dry run of the preparation, read from stdin, must plan no
rule besides the allowed ones.
--write-manifest: list the files under --root, with sizes, after preparing.
"""

import argparse
import os
import re
import subprocess
import sys

NOTHING_TO_DO = "Nothing to be done"
JOB_ROW = re.compile(r"^(\w+)\s+(\d+)$")


def planned_rules(dry_run_output):
    """Rule name to job count from the job table of a dry run; None without a plan."""
    if NOTHING_TO_DO in dry_run_output:
        return {}
    lines = dry_run_output.splitlines()
    if "Job stats:" not in lines:
        return None
    rules = {}
    for line in lines[lines.index("Job stats:") + 1:]:
        row = JOB_ROW.match(line.strip())
        if row and row.group(1) == "total":
            break
        if row:
            rules[row.group(1)] = int(row.group(2))
    return rules


def plan_problem(dry_run_output, allowed):
    rules = planned_rules(dry_run_output)
    if rules is None:
        return "the dry run printed no plan"
    redone = {rule: count for rule, count in rules.items() if rule not in allowed}
    if redone:
        listing = ", ".join(f"{rule} ({count})" for rule, count in sorted(redone.items()))
        return f"preparation would be redone: {listing}"
    return None


def commit_problem(commit_file, head):
    try:
        with open(commit_file) as handle:
            prepared = handle.read().strip()
    except FileNotFoundError:
        return f"{commit_file} is missing: the tree was not prepared with make euler"
    if prepared != head:
        return f"prepared at {prepared}, checked out at {head}"
    return None


def file_sizes(root):
    """Path relative to root to size in bytes, links included, not followed."""
    sizes = {}
    for directory, _, names in os.walk(root):
        for name in names:
            path = os.path.join(directory, name)
            sizes[os.path.relpath(path, root)] = os.lstat(path).st_size
    return sizes


def write_manifest(root, manifest):
    with open(manifest, "w") as handle:
        for path, size in sorted(file_sizes(root).items()):
            handle.write(f"{path}\t{size}\n")


def manifest_problem(root, manifest):
    try:
        with open(manifest) as handle:
            expected = dict(line.rstrip("\n").split("\t") for line in handle)
    except FileNotFoundError:
        return f"{manifest} is missing: the tree was not prepared with make euler"
    found = file_sizes(root)
    missing = sorted(path for path in expected if path not in found)
    resized = sorted(path for path in expected if path in found and found[path] != int(expected[path]))
    if missing or resized:
        examples = ", ".join((missing + resized)[:3])
        return f"{len(missing)} files missing and {len(resized)} of another size under {root}, e.g. {examples}"
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--commit-file")
    parser.add_argument("--root", help="directory the manifest describes")
    parser.add_argument("--manifest")
    parser.add_argument("--write-manifest", action="store_true")
    parser.add_argument("--plan", action="store_true")
    parser.add_argument("--allow", nargs="*", default=["build_fastder"])
    args = parser.parse_args()

    if args.write_manifest:
        write_manifest(args.root, args.manifest)
        return 0

    problems = []
    if args.manifest:
        problems.append(manifest_problem(args.root, args.manifest))
    if args.commit_file:
        head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True,
                              text=True, check=True).stdout.strip()
        problems.append(commit_problem(args.commit_file, head))
    if args.plan:
        problems.append(plan_problem(sys.stdin.read(), set(args.allow)))
    problems = [problem for problem in problems if problem]
    for problem in problems:
        print(f"check_prepared: {problem}", file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
