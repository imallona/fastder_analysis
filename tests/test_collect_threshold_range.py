"""Tests for collect_threshold_range.py.

Covers a locus called in the case group only, in both groups, in neither, a
separating range with a gap, and runs that differ in another parameter.
"""

import csv

from collect_threshold_range import called, collect, loci_per_threshold, read_loci, summarise

LOCUS = {"gene": "STMN2", "chrom": "chr8", "start": 1000, "end": 1200}


def write_gtf(group_dir, param_id, exons):
    run_dir = group_dir / param_id
    run_dir.mkdir(parents=True)
    lines = ["#format: gtf\n"]
    for chrom, start, end in exons:
        lines.append("\t".join([chrom, "fastder", "exon", str(start), str(end),
                                "1.0", ".", ".", 'gene_id "g"; transcript_id "t";']) + "\n")
    (run_dir / "output.gtf").write_text("".join(lines))


def make_groups(tmp_path, runs):
    """runs maps param_id to (case exons, control exons)."""
    case, control = tmp_path / "knockdown", tmp_path / "control"
    for param_id, (case_exons, control_exons) in runs.items():
        write_gtf(case, param_id, case_exons)
        write_gtf(control, param_id, control_exons)
    return str(case), str(control)


def test_overlap_is_inclusive_at_both_ends():
    assert called(LOCUS, {"chr8": [(1200, 1300)]})
    assert called(LOCUS, {"chr8": [(900, 1000)]})
    assert not called(LOCUS, {"chr8": [(1201, 1300)]})
    assert not called(LOCUS, {"chr19": [(1000, 1200)]})


def test_loci_use_the_cryptic_exon_columns(tmp_path):
    tsv = tmp_path / "loci.tsv"
    tsv.write_text("gene\tchrom\tstart\tend\tce_start\tce_end\tsource\n"
                   "STMN2\tchr8\t1\t5000\t1000\t1200\tsomewhere\n")
    assert read_loci(tsv) == [LOCUS]


def test_a_locus_called_in_both_groups_is_not_separated(tmp_path):
    hit = [("chr8", 1100, 1150)]
    case, control = make_groups(tmp_path, {
        "mc0.01": (hit, hit),
        "mc0.5": (hit, []),
        "mc5.0": ([], []),
    })
    rows = collect([LOCUS], case, control)
    assert {r["param_id"]: r["separates"] for r in rows} == {"mc0.01": 0, "mc0.5": 1, "mc5.0": 0}
    summary = summarise(rows)[0]
    assert (summary["lowest_separating"], summary["highest_separating"]) == (0.5, 0.5)
    assert summary["thresholds_tested"] == 3


def test_a_gap_in_the_range_is_reported(tmp_path):
    hit = [("chr8", 1100, 1150)]
    case, control = make_groups(tmp_path, {
        "mc0.1": (hit, []),
        "mc0.2": (hit, hit),
        "mc0.5": (hit, []),
    })
    summary = summarise(collect([LOCUS], case, control))[0]
    assert (summary["lowest_separating"], summary["highest_separating"]) == (0.1, 0.5)
    assert summary["contiguous"] == 0


def test_a_locus_never_separated_has_no_range(tmp_path):
    case, control = make_groups(tmp_path, {"mc0.1": ([], []), "mc0.5": ([], [])})
    summary = summarise(collect([LOCUS], case, control))[0]
    assert summary["thresholds_separating"] == 0
    assert summary["lowest_separating"] == ""


def test_other_parameters_keep_their_own_range(tmp_path):
    hit = [("chr8", 1100, 1150)]
    case, control = make_groups(tmp_path, {
        "mc0.1_ct1.0": (hit, []),
        "mc0.1_ct1000.0": (hit, hit),
    })
    summary = {s["other_params"]: s["thresholds_separating"]
               for s in summarise(collect([LOCUS], case, control))}
    assert summary == {"coverage_tolerance=1.0": 1, "coverage_tolerance=1000.0": 0}


def test_the_threshold_separating_most_loci_is_flagged(tmp_path):
    other = {"gene": "HDGFL2", "chrom": "chr19", "start": 500, "end": 600}
    stmn2, hdgfl2 = [("chr8", 1100, 1150)], [("chr19", 520, 560)]
    case, control = make_groups(tmp_path, {
        "mc0.002": (stmn2 + hdgfl2, stmn2),
        "mc0.005": (stmn2 + hdgfl2, []),
        "mc0.02": (stmn2, []),
        "mc0.5": ([], []),
    })
    counts = {row["min_coverage"]: row
              for row in loci_per_threshold(collect([LOCUS, other], case, control))}
    assert [counts[t]["loci_separating"] for t in (0.002, 0.005, 0.02, 0.5)] == [1, 2, 1, 0]
    assert counts[0.005]["loci"] == "HDGFL2 STMN2"
    assert [t for t, row in counts.items() if row["most_loci"]] == [0.005]


def test_the_shipped_loci_table_parses():
    from pathlib import Path
    tsv = Path(__file__).resolve().parents[1] / "config" / "tdp43_cryptic_exons.tsv"
    loci = read_loci(tsv)
    assert [locus["gene"] for locus in loci] == ["STMN2", "HDGFL2", "ELAVL3", "CELF5", "KCNQ2"]
    with open(tsv, newline="") as fh:
        assert all(row["source"] for row in csv.DictReader(fh, delimiter="\t"))
