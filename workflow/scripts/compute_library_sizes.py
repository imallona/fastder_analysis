#!/usr/bin/env python3
"""Per-sample library sizes for every BigWig in a directory.

The library size is the sum of length times value over the whole file, which
is what fastder accumulates in Parser.cpp and what recount3 calls the AUC. It
covers every chromosome present, not only the ones a run analyses, so a given
CPM threshold is the same absolute cutoff whatever subset is analysed.

Computing it here, once, rather than inside each tool keeps the runtime
comparison honest. fastder reads the total from the BigWig summary header for
free, while derfinder would have to import every chromosome to reach the same
number; charging that to derfinder inside its benchmarked rule would inflate
its wall time by the ratio of the genome to the analysed subset.

Output is a TSV of bigwig path, sample id and library size, sorted by path.
"""

import argparse
import os
import sys

import pyBigWig


def library_size(bw_path):
    """Sum of length times value over the whole BigWig.

    The summary header holds this sum already. A file written without one is
    walked chromosome by chromosome.
    """
    bw = pyBigWig.open(bw_path)
    try:
        total = float(bw.header().get("sumData", 0.0))
        if total > 0:
            return total
        total = 0.0
        for chrom in bw.chroms():
            for start, end, value in bw.intervals(chrom) or []:
                total += (end - start) * float(value)
    finally:
        bw.close()
    return total


def sample_id(bw_path):
    """Strip the .bw suffix and any strand marker fastder also strips."""
    name = os.path.basename(bw_path)
    for suffix in (".all.bw", ".unique.bw", ".plus.bw", ".minus.bw", ".bw"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bigwig-dir", required=True,
                        help="Directory holding the per-sample BigWig files")
    parser.add_argument("--out", required=True, help="Output TSV path")
    args = parser.parse_args()

    paths = sorted(os.path.join(args.bigwig_dir, f)
                   for f in os.listdir(args.bigwig_dir) if f.endswith(".bw"))
    if not paths:
        sys.exit(f"No .bw files in {args.bigwig_dir}")

    with open(args.out, "w") as out:
        out.write("bigwig\tsample\tlibrary_size\n")
        for path in paths:
            size = library_size(path)
            if size <= 0:
                print(f"WARN: {path} has an empty library size", file=sys.stderr)
            out.write(f"{path}\t{sample_id(path)}\t{size:.6f}\n")


if __name__ == "__main__":
    main()
