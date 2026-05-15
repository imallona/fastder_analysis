"""Subset a recount3 study junction matrix to one sample group.

recount3 ships per-study junction files that are already in the Monorail
RR/MM format: sra.junctions.<study>.ALL.RR (junction coordinates, one row
per junction), sra.junctions.<study>.ALL.MM (MatrixMarket sparse matrix,
junctions by samples) and sra.junctions.<study>.ALL.ID (the rail_id of each
MM column, in column order).

This script keeps only the junctions on the requested chromosomes and only
the MM columns belonging to the requested samples, then renumbers junction
ids and sample ids to contiguous 1-based ranges. The output matches what
emit_lean_mm_rr.py produces for the monorail_light backend and what fastder's
lean RR/MM parser expects.

The requested samples are given as SRA run accessions. The MM columns are
keyed by recount3 rail_id, so the study metadata table is used to map
rail_id to run accession.

Outputs:
    <out_prefix>.RR  TSV with header, columns chromosome start end length
                     strand annotated left_motif right_motif left_annotated
                     right_annotated. The two annotation columns are written
                     as "." because fastder reads them but never uses them.
    <out_prefix>.MM  MatrixMarket coordinate matrix, n_junctions by n_samples.
                     fastder uses MM as a presence mask, so each kept
                     (junction, sample) entry is written with a count of 1.
    <out_prefix>.samples.tsv  rail_id, sample_id, study_id, one row per
                     sample in MM column order. Consumed by
                     create_bigwig_list.sh.

Usage:
    python subset_recount3_junctions.py \\
        --rr study.ALL.RR --mm study.ALL.MM --id study.ALL.ID \\
        --metadata study.recount_project.tsv --study SRP166282 \\
        --sample SRR8083867 --sample SRR8083868 \\
        --chromosomes chr8 chr19 \\
        --out-prefix /path/junctions.ALL
"""
import argparse
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--rr", required=True, help="recount3 ALL.RR file (decompressed)")
    p.add_argument("--mm", required=True, help="recount3 ALL.MM file (decompressed)")
    p.add_argument("--id", required=True, help="recount3 ALL.ID file (decompressed)")
    p.add_argument("--metadata", required=True,
                   help="recount3 recount_project metadata TSV (decompressed)")
    p.add_argument("--study", required=True, help="SRA study accession")
    p.add_argument("--sample", action="append", required=True,
                   help="SRA run accession to keep (repeat; sets MM column order)")
    p.add_argument("--chromosomes", nargs="+", required=True,
                   help="Chromosome whitelist; junctions on others are dropped")
    p.add_argument("--out-prefix", required=True,
                   help="Output path prefix; writes <prefix>.RR, <prefix>.MM, "
                        "<prefix>.samples.tsv")
    return p.parse_args()


def read_railid_to_run(metadata_path):
    """Map recount3 rail_id to SRA run accession from the project metadata."""
    railid_to_run = {}
    with open(metadata_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rail_col = header.index("rail_id")
        run_col = header.index("external_id")
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            f = line.split("\t")
            railid_to_run[f[rail_col]] = f[run_col]
    return railid_to_run


def read_id_column_order(id_path):
    """Return the rail_id of each MM column, in column order (column 1 first)."""
    rail_ids = []
    with open(id_path) as fh:
        first = fh.readline().rstrip("\n")
        # The file starts with a "rail_id" header line.
        if first and first != "rail_id":
            rail_ids.append(first)
        for line in fh:
            line = line.rstrip("\n")
            if line:
                rail_ids.append(line)
    return rail_ids


def main():
    args = parse_args()
    chrom_set = set(args.chromosomes)

    railid_to_run = read_railid_to_run(args.metadata)
    column_rail_ids = read_id_column_order(args.id)
    # recount3 MM columns are 1-based and follow ID-file order.
    run_to_recount3_col = {}
    for col_index, rail_id in enumerate(column_rail_ids, start=1):
        run = railid_to_run.get(rail_id)
        if run is not None:
            run_to_recount3_col[run] = col_index

    missing = [s for s in args.sample if s not in run_to_recount3_col]
    if missing:
        sys.exit(f"[subset_recount3_junctions] runs not found in study "
                 f"{args.study} junction matrix: {', '.join(missing)}")

    # Output column order follows the --sample argument order.
    recount3_col_to_out_col = {
        run_to_recount3_col[run]: out_col
        for out_col, run in enumerate(args.sample, start=1)
    }
    kept_recount3_cols = set(recount3_col_to_out_col)

    # First pass over RR: decide which junctions to keep and renumber them.
    # recount3 RR row N (after the header) corresponds to MM junction id N.
    out_prefix = Path(args.out_prefix)
    rr_out_path = Path(str(out_prefix) + ".RR")
    mm_out_path = Path(str(out_prefix) + ".MM")
    samples_out_path = Path(str(out_prefix) + ".samples.tsv")

    recount3_sjid_to_out_sjid = {}
    with open(args.rr) as rr_in, open(rr_out_path, "w") as rr_out:
        rr_out.write("chromosome\tstart\tend\tlength\tstrand\tannotated"
                     "\tleft_motif\tright_motif\tleft_annotated\tright_annotated\n")
        header = rr_in.readline()
        if not header.startswith("chromosome"):
            sys.exit("[subset_recount3_junctions] unexpected RR header: "
                     + header.rstrip("\n"))
        recount3_sjid = 0
        out_sjid = 0
        for line in rr_in:
            line = line.rstrip("\n")
            if not line:
                continue
            recount3_sjid += 1
            f = line.split("\t")
            if f[0] not in chrom_set:
                continue
            out_sjid += 1
            recount3_sjid_to_out_sjid[recount3_sjid] = out_sjid
            # Pass through coordinate and motif columns, suppress the two
            # annotation columns to keep the file lean.
            rr_out.write("\t".join(f[0:8]) + "\t.\t.\n")

    n_junctions = out_sjid

    # Second pass over MM: keep entries whose junction and column survived.
    triples = []
    with open(args.mm) as mm_in:
        dims_seen = False
        for line in mm_in:
            line = line.strip()
            if not line or line.startswith("%"):
                continue
            if not dims_seen:
                # The first non-comment line is "n_rows n_cols n_nonzero".
                dims_seen = True
                continue
            parts = line.split()
            recount3_sjid = int(parts[0])
            recount3_col = int(parts[1])
            if recount3_col not in kept_recount3_cols:
                continue
            out_sjid = recount3_sjid_to_out_sjid.get(recount3_sjid)
            if out_sjid is None:
                continue
            out_col = recount3_col_to_out_col[recount3_col]
            triples.append((out_sjid, out_col))

    n_samples = len(args.sample)
    with open(mm_out_path, "w") as mm_out:
        mm_out.write("%%MatrixMarket matrix coordinate integer general\n")
        mm_out.write("%-----------------------------------------------\n")
        mm_out.write(f"{n_junctions}\t{n_samples}\t{len(triples)}\n")
        for out_sjid, out_col in triples:
            mm_out.write(f"{out_sjid}\t{out_col}\t1\n")

    with open(samples_out_path, "w") as samples_out:
        samples_out.write("rail_id\tsample_id\tstudy_id\n")
        for out_col, run in enumerate(args.sample, start=1):
            samples_out.write(f"{out_col}\t{run}\t{args.study}\n")

    print(f"[subset_recount3_junctions] study {args.study}: "
          f"junctions {n_junctions}  samples {n_samples}  "
          f"nonzeros {len(triples)}", file=sys.stderr)


if __name__ == "__main__":
    main()
