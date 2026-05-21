"""Naive coverage segmenter as a baseline for fastder.

Reads every BigWig in --bigwig-dir, normalizes each per-sample coverage
track to CPM (counts per million) using the same library_size formula
fastder applies internally, averages the per-sample CPMs across all
samples (zero-included for samples without coverage at a position), and
emits one GTF transcript per maximal run of bases whose mean CPM is at
or above --cutoff.

This deliberately mirrors what fastder's Averager does, minus the SJ-
aware stitching step. The point is to isolate the lift fastder gets
from its SJ-aware stitching over a thresholded coverage segmenter:
identical input data, identical normalization, identical aggregation,
only the region-call step differs.

Parameter equivalence with fastder:
  --cutoff      <-> fastder --min-coverage   (CPM)
  --min-length  <-> fastder --min-length     (bp; post-filter)
  --chromosomes <-> fastder --chr            (used for library_size scope)
"""
import argparse
import glob
import os
import os.path as op

import pyBigWig
import numpy as np


def load_bigwig_samples(bw_dir):
    """Return one entry per sample, each a list of that sample's BigWig
    files. An unstranded sample has a single .all.bw. A stranded sample has
    a .plus.bw and a .minus.bw; those are two strand tracks of the same
    sample and must be treated as one sample, not two, or the CPM library
    sizes and the averaging run over twice as many pseudo-samples.

    TODO(deprecate): the stranded branch supports the workflow's
    stranded=true hack (see config.yaml). No paper figure uses it, but
    removing it needs the original stranded-path contributor's
    agreement. Until then we keep both branches.
    """
    all_bw = sorted(glob.glob(op.join(bw_dir, "*.all.bw")))
    if all_bw:
        return [[p] for p in all_bw]
    stranded = (glob.glob(op.join(bw_dir, "*.plus.bw"))
                + glob.glob(op.join(bw_dir, "*.minus.bw")))
    if not stranded:
        raise FileNotFoundError(f"no .bw files in {bw_dir}")
    by_sample = {}
    for path in stranded:
        name = op.basename(path)
        for suffix in (".plus.bw", ".minus.bw"):
            if name.endswith(suffix):
                by_sample.setdefault(name[: -len(suffix)], []).append(path)
                break
    return [sorted(by_sample[s]) for s in sorted(by_sample)]


def common_chroms(bw_paths):
    chroms = None
    for p in bw_paths:
        bw = pyBigWig.open(p)
        keys = set(bw.chroms().keys())
        chroms = keys if chroms is None else (chroms & keys)
        bw.close()
    return sorted(chroms or [])


def chrom_length(bw_paths, chrom):
    for p in bw_paths:
        bw = pyBigWig.open(p)
        try:
            length = bw.chroms().get(chrom)
        finally:
            bw.close()
        if length:
            return length
    raise RuntimeError(f"chromosome {chrom} not in any BigWig")


def library_size(bw_path, chroms):
    """Sum of length * value over the given chromosomes for one BigWig.

    Matches fastder's library_size accumulation in Parser.cpp: per-row
    `total_reads = (end - start) * coverage`, summed over only the
    chromosomes the user passed via --chr. With per-base coverage we get
    the same result by computing length * value for each per-base record
    coming out of pyBigWig.intervals().
    """
    bw = pyBigWig.open(bw_path)
    try:
        bw_chroms = set(bw.chroms().keys())
        total = 0.0
        for chrom in chroms:
            if chrom not in bw_chroms:
                continue
            intervals = bw.intervals(chrom)
            if not intervals:
                continue
            for start, end, value in intervals:
                total += (end - start) * float(value)
    finally:
        bw.close()
    return total


def mean_cpm_coverage(samples, chrom, length, cpm_factors):
    """Per-base mean CPM across samples.

    Each sample is a list of BigWig files: one for an unstranded sample,
    two (plus and minus) for a stranded sample. A sample's raw coverage is
    the sum of its tracks. cpm_factors[i] = library_size_i / 1e6 turns that
    into CPM. Per-sample CPM is then averaged across samples (a sample with
    no coverage at a position contributes zero, matching fastder's Averager).
    """
    if len(samples) != len(cpm_factors):
        raise ValueError("samples and cpm_factors must have the same length")
    acc = np.zeros(length, dtype=np.float64)
    for paths, factor in zip(samples, cpm_factors):
        if factor <= 0:
            continue
        sample_cov = np.zeros(length, dtype=np.float64)
        for path in paths:
            bw = pyBigWig.open(path)
            try:
                vals = bw.values(chrom, 0, length, numpy=True)
            finally:
                bw.close()
            if vals is None:
                continue
            np.nan_to_num(vals, copy=False, nan=0.0)
            sample_cov += vals
        acc += sample_cov / factor
    return acc / len(samples)


def find_regions(coverage, cutoff, min_length):
    """Maximal runs of consecutive bases with coverage at or above cutoff."""
    above = coverage >= cutoff
    if not above.any():
        return []
    flips = np.diff(above.astype(np.int8))
    starts = np.where(flips == 1)[0] + 1
    ends = np.where(flips == -1)[0] + 1
    if above[0]:
        starts = np.r_[0, starts]
    if above[-1]:
        ends = np.r_[ends, len(coverage)]
    regions = []
    for s, e in zip(starts, ends):
        if e - s < min_length:
            continue
        regions.append((int(s), int(e), float(coverage[s:e].mean())))
    return regions


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bigwig-dir", required=True)
    ap.add_argument("--out-gtf", required=True)
    ap.add_argument("--cutoff", type=float, default=0.05,
                    help="Coverage threshold in CPM (matches fastder --min-coverage)")
    ap.add_argument("--min-length", type=int, default=10)
    ap.add_argument("--chromosomes", nargs="+", default=None,
                    help="Chromosomes to analyse and to scope library_size to. "
                         "Default = intersection across all BigWigs.")
    args = ap.parse_args()

    samples = load_bigwig_samples(args.bigwig_dir)
    flat_paths = [p for sample in samples for p in sample]
    chroms = args.chromosomes if args.chromosomes else common_chroms(flat_paths)

    # CPM scaling factor per sample, restricted to the user's chromosome set.
    # A stranded sample's library size is the sum of its plus and minus tracks.
    cpm_factors = [sum(library_size(p, chroms) for p in sample) / 1e6
                   for sample in samples]
    for sample, factor in zip(samples, cpm_factors):
        if factor <= 0:
            print(f"[megadepth_baseline] WARN: {op.basename(sample[0])} has empty "
                  f"library_size on chromosomes {chroms}; sample will be skipped.")

    os.makedirs(op.dirname(args.out_gtf), exist_ok=True)
    with open(args.out_gtf, "w") as out:
        out.write("# megadepth_baseline: CPM-normalized mean-coverage segmenter\n")
        out.write(f"# cutoff={args.cutoff} CPM, min_length={args.min_length} bp, "
                  f"chromosomes={','.join(chroms)}\n")
        gid = 0
        for chrom in chroms:
            length = chrom_length(flat_paths, chrom)
            cov = mean_cpm_coverage(samples, chrom, length, cpm_factors)
            regions = find_regions(cov, args.cutoff, args.min_length)
            for start, end, mean_cov in regions:
                gid += 1
                gene = f"gene{gid}"
                tx = f"tx{gid}"
                gtf_start = start + 1
                gtf_end = end
                attrs_gene = f'gene_id "{gene}"; gene_name "{gene}_naive";'
                attrs_tx = f'gene_id "{gene}"; transcript_id "{tx}";'
                attrs_exon = f'gene_id "{gene}"; transcript_id "{tx}"; exon_number "1";'
                out.write(f'{chrom}\tmegadepth_baseline\tgene\t{gtf_start}\t{gtf_end}\t{mean_cov:.4f}\t.\t.\t{attrs_gene}\n')
                out.write(f'{chrom}\tmegadepth_baseline\ttranscript\t{gtf_start}\t{gtf_end}\t{mean_cov:.4f}\t.\t.\t{attrs_tx}\n')
                out.write(f'{chrom}\tmegadepth_baseline\texon\t{gtf_start}\t{gtf_end}\t{mean_cov:.4f}\t.\t.\t{attrs_exon}\n')


if __name__ == "__main__":
    main()
