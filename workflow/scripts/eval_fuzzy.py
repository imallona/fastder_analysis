"""Fuzzy fastder vs reference comparison without gffcompare's exact-boundary requirement. Produces three CSVs: reciprocal-best Jaccard per ref transcript, per-fastder-exon distance to the nearest ref exon (same strand), and locus-recall at fixed coverage thresholds. All bedtools calls use -s for strand-aware overlap."""
import argparse
import csv
import os.path as op
import re
import subprocess
import sys
import tempfile
from collections import defaultdict


GTF_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')
GFF_ATTR_RE = re.compile(r'(\w+)=([^;]+)')


def parse_gtf_attributes(attr_str):
    return dict(GTF_ATTR_RE.findall(attr_str))


def parse_gff_attributes(attr_str):
    return dict(GFF_ATTR_RE.findall(attr_str))


def gtf_or_gff_to_bed(in_path, out_path, attribute_parser):
    """Write BED-6: chrom, start (0-based), end, transcript_id, gene_id, strand."""
    with open(in_path) as src, open(out_path, "w") as dst:
        for line in src:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "exon":
                continue
            chrom, start, end, strand = cols[0], int(cols[3]) - 1, int(cols[4]), cols[6]
            attrs = attribute_parser(cols[8])
            tx = attrs.get("transcript_id", "")
            gene = attrs.get("gene_id", "")
            if not tx or not gene:
                continue
            dst.write(f"{chrom}\t{start}\t{end}\t{tx}\t{gene}\t{strand}\n")


def total_bp_per_key(bed_path, key_col):
    """Sum (end - start) grouped by the column at key_col (0-indexed)."""
    totals = defaultdict(int)
    with open(bed_path) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            totals[cols[key_col]] += int(cols[2]) - int(cols[1])
    return totals


def reciprocal_best_jaccard(ref_bed, fastder_bed, sample, param_id, scenario, out_path):
    ref_total = total_bp_per_key(ref_bed, 3)
    fastder_total = total_bp_per_key(fastder_bed, 3)
    if not ref_total:
        with open(out_path, "w") as f:
            csv.writer(f).writerow([
                "scenario", "sample", "param_id", "ref_transcript", "ref_gene",
                "fastder_transcript", "overlap_bp", "ref_bp", "fastder_bp", "jaccard"
            ])
        return

    pair_overlap = defaultdict(int)
    pair_meta = {}
    cmd = ["bedtools", "intersect", "-s", "-wo", "-a", ref_bed, "-b", fastder_bed]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    for line in proc.stdout.splitlines():
        cols = line.split("\t")
        ref_tx = cols[3]
        ref_gene = cols[4]
        fastder_tx = cols[9]
        overlap = int(cols[12])
        pair_overlap[(ref_tx, fastder_tx)] += overlap
        pair_meta[(ref_tx, fastder_tx)] = ref_gene

    best = {}
    for (ref_tx, fastder_tx), overlap in pair_overlap.items():
        union = ref_total[ref_tx] + fastder_total.get(fastder_tx, 0) - overlap
        if union <= 0:
            continue
        jaccard = overlap / union
        cur = best.get(ref_tx)
        if cur is None or jaccard > cur["jaccard"]:
            best[ref_tx] = {
                "fastder_transcript": fastder_tx,
                "overlap_bp": overlap,
                "fastder_bp": fastder_total.get(fastder_tx, 0),
                "jaccard": jaccard,
                "ref_gene": pair_meta[(ref_tx, fastder_tx)],
            }

    with open(out_path, "w") as f:
        w = csv.writer(f)
        w.writerow([
            "scenario", "sample", "param_id", "ref_transcript", "ref_gene",
            "fastder_transcript", "overlap_bp", "ref_bp", "fastder_bp", "jaccard"
        ])
        for ref_tx, ref_bp in ref_total.items():
            entry = best.get(ref_tx)
            if entry is None:
                w.writerow([scenario, sample, param_id, ref_tx, "", "", 0, ref_bp, 0, 0.0])
            else:
                w.writerow([
                    scenario, sample, param_id, ref_tx, entry["ref_gene"],
                    entry["fastder_transcript"], entry["overlap_bp"], ref_bp,
                    entry["fastder_bp"], f"{entry['jaccard']:.6f}",
                ])


def boundary_distances(ref_bed, fastder_bed, sample, param_id, scenario, out_path):
    """For each fastder exon, signed distance from its start to nearest ref
    exon start (same strand) and from its end to nearest ref exon end. The
    closest call is run twice on single-position BEDs so each end is matched
    independently. Negative distance = fastder boundary is upstream of the
    nearest ref boundary on the feature's own strand.
    """
    def project(bed_in, bed_out, which):
        with open(bed_in) as src, open(bed_out, "w") as dst:
            for line in src:
                cols = line.rstrip("\n").split("\t")
                start = int(cols[1])
                end = int(cols[2])
                pos = start if which == "start" else end - 1
                dst.write(f"{cols[0]}\t{pos}\t{pos + 1}\t{cols[3]}\t{cols[4]}\t{cols[5]}\n")

    rows = []
    with tempfile.TemporaryDirectory() as td:
        for which in ("start", "end"):
            ref_one = op.join(td, f"ref_{which}.bed")
            fast_one = op.join(td, f"fast_{which}.bed")
            project(ref_bed, ref_one, which)
            project(fastder_bed, fast_one, which)
            for path in (ref_one, fast_one):
                subprocess.run(["sort", "-k1,1", "-k2,2n", path, "-o", path], check=True)
            cmd = ["bedtools", "closest", "-s", "-D", "ref", "-t", "first",
                   "-a", fast_one, "-b", ref_one]
            proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
            for line in proc.stdout.splitlines():
                cols = line.split("\t")
                if len(cols) < 13 or cols[6] == ".":
                    continue
                fast_chrom = cols[0]
                fast_pos = int(cols[1])
                fast_tx = cols[3]
                fast_strand = cols[5]
                ref_tx = cols[9]
                distance = int(cols[12])
                rows.append({
                    "scenario": scenario,
                    "sample": sample,
                    "param_id": param_id,
                    "boundary": which,
                    "chrom": fast_chrom,
                    "fastder_transcript": fast_tx,
                    "fastder_pos": fast_pos,
                    "ref_transcript": ref_tx,
                    "distance": distance,
                    "strand": fast_strand,
                })

    with open(out_path, "w") as f:
        if rows:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            for r in rows:
                w.writerow(r)
        else:
            f.write("scenario,sample,param_id,boundary,chrom,fastder_transcript,"
                    "fastder_pos,ref_transcript,distance,strand\n")


def locus_recall(ref_bed, fastder_bed, sample, param_id, scenario, out_path,
                 thresholds=tuple(round(0.05 * i, 2) for i in range(1, 21))):
    """Per ref locus, what fraction of its exonic bp is covered by any
    fastder exon on the same strand. Recall at threshold f = fraction of
    loci with covered/total >= f.
    """
    locus_total = total_bp_per_key(ref_bed, 4)
    if not locus_total:
        with open(out_path, "w") as f:
            csv.writer(f).writerow(
                ["scenario", "sample", "param_id", "threshold", "n_loci_recovered", "n_loci_total", "recall"])
        return

    cmd = ["bedtools", "intersect", "-s", "-wo", "-a", ref_bed, "-b", fastder_bed]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    locus_covered_intervals = defaultdict(list)
    for line in proc.stdout.splitlines():
        cols = line.split("\t")
        gene = cols[4]
        ref_start = int(cols[1])
        ref_end = int(cols[2])
        f_start = int(cols[7])
        f_end = int(cols[8])
        s = max(ref_start, f_start)
        e = min(ref_end, f_end)
        if e > s:
            locus_covered_intervals[gene].append((s, e))

    locus_covered = {}
    for gene, intervals in locus_covered_intervals.items():
        intervals.sort()
        merged_bp = 0
        cur_s, cur_e = intervals[0]
        for s, e in intervals[1:]:
            if s <= cur_e:
                if e > cur_e:
                    cur_e = e
            else:
                merged_bp += cur_e - cur_s
                cur_s, cur_e = s, e
        merged_bp += cur_e - cur_s
        locus_covered[gene] = merged_bp

    n_total = len(locus_total)
    with open(out_path, "w") as f:
        w = csv.writer(f)
        w.writerow(["scenario", "sample", "param_id", "threshold", "n_loci_recovered", "n_loci_total", "recall"])
        for thr in thresholds:
            n_rec = sum(
                1 for g, total in locus_total.items()
                if total > 0 and locus_covered.get(g, 0) / total >= thr
            )
            w.writerow([scenario, sample, param_id, thr, n_rec, n_total, f"{n_rec / n_total:.6f}"])


def strand_concordance(ref_bed, fastder_bed, sample, param_id, scenario, out_path):
    """Per fastder transcript, find the best-overlapping reference transcript
    on any strand and classify the strand outcome:
      concordant: fastder strand matches the best ref partner's strand
      discordant: fastder strand opposite of best ref partner's strand
      unstranded: fastder transcript carries strand '.'
      unmatched: fastder transcript has no ref overlap at all
    """
    fastder_strand = {}
    with open(fastder_bed) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            fastder_strand[cols[3]] = cols[5]

    cmd = ["bedtools", "intersect", "-wo", "-a", fastder_bed, "-b", ref_bed]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)

    pair_overlap = defaultdict(int)
    for line in proc.stdout.splitlines():
        cols = line.split("\t")
        fastder_tx = cols[3]
        ref_strand = cols[11]
        overlap = int(cols[12])
        pair_overlap[(fastder_tx, ref_strand)] += overlap

    best_ref_strand = {}
    for (fastder_tx, ref_strand), overlap in pair_overlap.items():
        cur = best_ref_strand.get(fastder_tx)
        if cur is None or overlap > cur[1]:
            best_ref_strand[fastder_tx] = (ref_strand, overlap)

    counts = {"concordant": 0, "discordant": 0, "unstranded": 0, "unmatched": 0}
    for tx, fstrand in fastder_strand.items():
        if fstrand == ".":
            counts["unstranded"] += 1
        elif tx not in best_ref_strand:
            counts["unmatched"] += 1
        else:
            ref_strand = best_ref_strand[tx][0]
            if ref_strand == fstrand:
                counts["concordant"] += 1
            else:
                counts["discordant"] += 1

    with open(out_path, "w") as f:
        w = csv.writer(f)
        w.writerow(["scenario", "sample", "param_id", "category", "n_fastder_transcripts"])
        for cat in ("concordant", "discordant", "unstranded", "unmatched"):
            w.writerow([scenario, sample, param_id, cat, counts[cat]])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-gff", required=True)
    ap.add_argument("--fastder-gtf", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--param-id", required=True)
    ap.add_argument("--scenario", required=True)
    ap.add_argument("--out-jaccard", required=True)
    ap.add_argument("--out-distances", required=True)
    ap.add_argument("--out-recall", required=True)
    ap.add_argument("--out-strand", required=True)
    args = ap.parse_args()

    ref_is_gff3 = args.ref_gff.lower().endswith(".gff3") or args.ref_gff.lower().endswith(".gff")
    ref_parser = parse_gff_attributes if ref_is_gff3 else parse_gtf_attributes

    with tempfile.TemporaryDirectory() as td:
        ref_bed = op.join(td, "ref.bed")
        fastder_bed = op.join(td, "fastder.bed")
        gtf_or_gff_to_bed(args.ref_gff, ref_bed, ref_parser)
        gtf_or_gff_to_bed(args.fastder_gtf, fastder_bed, parse_gtf_attributes)
        for path in (ref_bed, fastder_bed):
            subprocess.run(["sort", "-k1,1", "-k2,2n", path, "-o", path], check=True)

        reciprocal_best_jaccard(ref_bed, fastder_bed, args.sample, args.param_id, args.scenario, args.out_jaccard)
        boundary_distances(ref_bed, fastder_bed, args.sample, args.param_id, args.scenario, args.out_distances)
        locus_recall(ref_bed, fastder_bed, args.sample, args.param_id, args.scenario, args.out_recall)
        strand_concordance(ref_bed, fastder_bed, args.sample, args.param_id, args.scenario, args.out_strand)


if __name__ == "__main__":
    main()
