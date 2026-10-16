"""Write each exon of a GTF as its own single-exon record.

An exon keeps its coordinates, score and strand. Its gene and transcript ids
are the original ones with the exon's position in the record appended.

Usage:
    python split_chains.py --gtf <in.gtf> --out <out.gtf>
"""
import argparse
import re
from collections import defaultdict

ATTRIBUTE = re.compile(r'(\w+) "([^"]*)"')


def split_records(lines):
    """Gene, transcript and exon lines for every exon line, in input order."""
    seen = defaultdict(int)
    for line in lines:
        if line.startswith("#") or not line.strip():
            yield line
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9 or cols[2] != "exon":
            continue
        attributes = dict(ATTRIBUTE.findall(cols[8]))
        transcript = attributes["transcript_id"]
        seen[transcript] += 1
        suffix = seen[transcript]
        gene_id = f'{attributes.get("gene_id", transcript)}.{suffix}'
        transcript_id = f"{transcript}.{suffix}"
        for feature, attribute_text in (
            ("gene", f'gene_id "{gene_id}";'),
            ("transcript", f'gene_id "{gene_id}"; transcript_id "{transcript_id}";'),
            ("exon", f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; exon_number "1";'),
        ):
            yield "\t".join([*cols[:2], feature, *cols[3:8], attribute_text]) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    with open(args.gtf) as src, open(args.out, "w") as out:
        out.writelines(split_records(src))


if __name__ == "__main__":
    main()
