"""IMPORTANT: produces per-scenario asimulator outputs, where scenario controls whether template + variant transcripts both contribute reads (template_and_variant, the ASimulatoR default) or only the variant transcripts do (variant_only, which removes both the template reads from the FASTQ and the template entries from the GFF so the truth set carries only the alternative isoforms)."""
import argparse
import gzip
import os
import os.path as op
import re
import shutil
import sys


GFF_ATTR_RE = re.compile(r'(\w+)=([^;]+)')


def parse_gff_attrs(attr_str):
    return dict(GFF_ATTR_RE.findall(attr_str))


def template_transcript_ids(gff_path):
    """Transcript IDs marked template=TRUE in the splicing_variants GFF."""
    template_ids = set()
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "transcript":
                continue
            attrs = parse_gff_attrs(cols[8])
            if attrs.get("template", "").upper() == "TRUE":
                tx = attrs.get("transcript_id", "")
                if tx:
                    template_ids.add(tx)
    return template_ids


def filter_gff(gff_in, gff_out, template_ids):
    """Drop transcript and exon rows whose transcript_id is a template. Drop
    gene rows that no longer have any non-template child. Lines without
    transcript_id (e.g. gene rows themselves) are handled in a second pass:
    keep a gene only if at least one of its variant transcripts survived."""
    surviving_gene_ids = set()
    surviving_lines = []
    with open(gff_in) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                surviving_lines.append((None, line))
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                surviving_lines.append((None, line))
                continue
            attrs = parse_gff_attrs(cols[8])
            tx = attrs.get("transcript_id", "")
            gene = attrs.get("gene_id", "")
            ftype = cols[2]
            if ftype == "gene":
                surviving_lines.append(("gene", line))
                continue
            if tx in template_ids:
                continue
            surviving_lines.append((ftype, line))
            if gene:
                surviving_gene_ids.add(gene)

    with open(gff_out, "w") as out:
        for tag, line in surviving_lines:
            if tag == "gene":
                cols = line.rstrip("\n").split("\t")
                attrs = parse_gff_attrs(cols[8])
                if attrs.get("gene_id", "") in surviving_gene_ids:
                    out.write(line)
            else:
                out.write(line)


READ_HEADER_TX_RE = re.compile(r'^@[^/]+/([^;\s]+)')


def fastq_iter(path):
    """Yield (header, seq, plus, qual) records. Handles plain or gz inputs."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        while True:
            h = f.readline()
            if not h:
                return
            s = f.readline()
            p = f.readline()
            q = f.readline()
            if not q:
                return
            yield h, s, p, q


def filter_fastq(fq_in, fq_out, template_ids):
    """Polyester (used by ASimulatoR) puts the originating transcript id in
    the read header, prefixed with @ and ending at the first / or :. Drop
    records whose transcript_id is a template. Output is uncompressed FASTQ
    to match ASimulatoR's default output format."""
    written = 0
    skipped = 0
    with open(fq_out, "w") as out:
        for h, s, p, q in fastq_iter(fq_in):
            m = READ_HEADER_TX_RE.match(h)
            tx = m.group(1) if m else ""
            if tx in template_ids:
                skipped += 1
                continue
            out.write(h)
            out.write(s)
            out.write(p)
            out.write(q)
            written += 1
    print(f"[make_scenario] {fq_in}: wrote {written}, dropped {skipped}", file=sys.stderr)


def passthrough(src, dst):
    if op.isfile(dst) or op.islink(dst):
        os.remove(dst)
    os.makedirs(op.dirname(dst), exist_ok=True)
    os.symlink(op.abspath(src), dst)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scenario", required=True,
                    choices=["template_and_variant", "variant_only"])
    ap.add_argument("--gff-in", required=True)
    ap.add_argument("--fq1-in", required=True)
    ap.add_argument("--fq2-in", required=True)
    ap.add_argument("--gff-out", required=True)
    ap.add_argument("--fq1-out", required=True)
    ap.add_argument("--fq2-out", required=True)
    args = ap.parse_args()

    os.makedirs(op.dirname(args.gff_out), exist_ok=True)

    if args.scenario == "template_and_variant":
        passthrough(args.gff_in, args.gff_out)
        passthrough(args.fq1_in, args.fq1_out)
        passthrough(args.fq2_in, args.fq2_out)
        return

    template_ids = template_transcript_ids(args.gff_in)
    print(f"[make_scenario] {len(template_ids)} template transcripts to drop",
          file=sys.stderr)
    filter_gff(args.gff_in, args.gff_out, template_ids)
    filter_fastq(args.fq1_in, args.fq1_out, template_ids)
    filter_fastq(args.fq2_in, args.fq2_out, template_ids)


if __name__ == "__main__":
    main()
