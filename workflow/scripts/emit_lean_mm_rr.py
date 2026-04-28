"""Aggregate STAR per-sample SJ.out.tab files into lean MM/RR files for fastder.

Outputs:
    <out_prefix>.RR  TSV with columns: chrom, start, end, length, strand,
                     annotated, left_motif, right_motif, left_annotated,
                     right_annotated. Annotation columns are emitted as "."
                     because fastder reads them into SJRow but never uses them
                     (verified against Parser.cpp / Integrator.cpp / main.cpp /
                     StitchedER.cpp). Suppressing the recount3 annotation
                     strings here drops ~200 bytes per junction off fastder's
                     resident set.

    <out_prefix>.MM  MatrixMarket coordinate-format sparse matrix
                     (n_junctions x n_samples).  fastder uses MM only as a
                     presence mask (Parser.cpp:213 stores sj_id but ignores
                     count), so a count of 1 per (sj, sample) pair is enough.

Junctions are filtered to the chromosomes listed in --chromosomes (matching
fastder's --chr restriction).  Restricting RR to the analysed chromosomes
shrinks rr_all_sj from ~9.5M (full hg38) to ~3.5k (chr21), since fastder
indexes that vector by sj_id and we re-number MM sj_ids to align.

Usage:
    python emit_lean_mm_rr.py \\
        --sample NAME --sj NAME.SJ.out.tab \\
        [--sample N2 --sj N2.SJ.out.tab ...] \\
        --chromosomes chr21 chr22 \\
        --min-reads 1 \\
        --out-prefix /path/junctions.ALL
"""
import argparse
import sys
from pathlib import Path


# STAR encodes the splice motif as an integer in column 5 of SJ.out.tab.
# Map it to the (left_dinuc, right_dinuc) pair stored in monorail's RR file.
STAR_MOTIF = {
    0: (".", "."),   # non-canonical
    1: ("GT", "AG"),
    2: ("CT", "AC"),
    3: ("GC", "AG"),
    4: ("CT", "GC"),
    5: ("AT", "AC"),
    6: ("GT", "AT"),
}
STAR_STRAND = {0: ".", 1: "+", 2: "-"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sample", action="append", required=True,
                   help="Sample name (repeat for each sample, in MM column order)")
    p.add_argument("--sj", action="append", required=True,
                   help="STAR SJ.out.tab path (repeat, must match --sample order)")
    p.add_argument("--chromosomes", nargs="+", required=True,
                   help="Chromosome whitelist; junctions on others are dropped")
    p.add_argument("--min-reads", type=int, default=1,
                   help="Skip a (sample, junction) pair when unique+multi reads < this")
    p.add_argument("--out-prefix", required=True,
                   help="Output path prefix; writes <prefix>.RR and <prefix>.MM")
    args = p.parse_args()
    if len(args.sample) != len(args.sj):
        p.error("--sample and --sj must be paired")
    return args


def read_sj_tab(path, chrom_set):
    """Yield (chrom, start, end, strand_code, motif_code, annotated, reads) tuples."""
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            chrom = f[0]
            if chrom not in chrom_set:
                continue
            start = int(f[1])
            end = int(f[2])
            strand_code = int(f[3])
            motif_code = int(f[4])
            annotated = int(f[5])
            reads = int(f[6]) + int(f[7])  # unique + multi
            yield (chrom, start, end, strand_code, motif_code, annotated, reads)


def main():
    args = parse_args()
    chrom_set = set(args.chromosomes)
    samples = args.sample
    sjs = args.sj

    # First pass: union of junction keys across samples, plus per-key metadata
    # from the first sample that reports it. Junction key = (chrom,start,end,strand).
    junctions = {}
    sample_calls = []  # list of dict[key -> reads], one per sample, parallel to samples

    for sj_path in sjs:
        calls = {}
        for chrom, start, end, strand_code, motif_code, annotated, reads in read_sj_tab(sj_path, chrom_set):
            key = (chrom, start, end, strand_code)
            if key not in junctions:
                junctions[key] = (motif_code, annotated)
            calls[key] = reads
        sample_calls.append(calls)

    # Sort junctions deterministically; sj_id is 1-based to match fastder's expectation
    sorted_keys = sorted(junctions.keys(), key=lambda k: (k[0], k[1], k[2], k[3]))
    key_to_sjid = {k: i + 1 for i, k in enumerate(sorted_keys)}

    out_prefix = Path(args.out_prefix)
    rr_path = out_prefix.with_suffix(out_prefix.suffix + ".RR") if out_prefix.suffix else Path(str(out_prefix) + ".RR")
    mm_path = out_prefix.with_suffix(out_prefix.suffix + ".MM") if out_prefix.suffix else Path(str(out_prefix) + ".MM")

    # Write RR (TSV with header). Columns must match fastder's SJRow operator>>
    # parser: chrom start end length strand annotated left_motif right_motif
    # left_annotated right_annotated. We emit "." for the unused annotation cols.
    with open(rr_path, "w") as out:
        out.write("chromosome\tstart\tend\tlength\tstrand\tannotated"
                  "\tleft_motif\tright_motif\tleft_annotated\tright_annotated\n")
        for key in sorted_keys:
            chrom, start, end, strand_code = key
            motif_code, annotated = junctions[key]
            length = end - start + 1
            strand = STAR_STRAND.get(strand_code, ".")
            left_motif, right_motif = STAR_MOTIF.get(motif_code, (".", "."))
            out.write(
                f"{chrom}\t{start}\t{end}\t{length}\t{strand}\t{annotated}"
                f"\t{left_motif}\t{right_motif}\t.\t.\n"
            )

    # Build MM triples (sj_id, sample_id, count) with sj_id 1-based and
    # sample_id 1-based. Apply the min-reads threshold per (sample, junction).
    triples = []
    for sample_idx, calls in enumerate(sample_calls, start=1):
        for key, reads in calls.items():
            if reads < args.min_reads:
                continue
            sj_id = key_to_sjid[key]
            # fastder ignores the count value — write 1 to keep MM compact
            triples.append((sj_id, sample_idx, 1))

    # MatrixMarket header: n_rows n_cols n_nonzero
    n_rows = len(sorted_keys)
    n_cols = len(samples)
    n_nz = len(triples)

    with open(mm_path, "w") as out:
        out.write("%%MatrixMarket matrix coordinate integer general\n")
        out.write("%-----------------------------------------------\n")
        out.write(f"{n_rows}\t{n_cols}\t{n_nz}\n")
        # MatrixMarket convention is column-major, but fastder treats MM as a
        # sparse presence list and doesn't require sorted order.
        for sj_id, sample_id, count in triples:
            out.write(f"{sj_id}\t{sample_id}\t{count}\n")

    print(f"[emit_lean_mm_rr] junctions: {n_rows}  samples: {n_cols}  nonzeros: {n_nz}",
          file=sys.stderr)
    print(f"[emit_lean_mm_rr] RR: {rr_path}", file=sys.stderr)
    print(f"[emit_lean_mm_rr] MM: {mm_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
