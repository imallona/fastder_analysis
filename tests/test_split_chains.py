"""Tests for split_chains.py.

Covers a three-exon chain, a single-exon record, and header lines.
"""

from split_chains import split_records


def gtf_line(feature, start, end, attributes, strand="+"):
    return "\t".join(["chr21", "fastder", feature, str(start), str(end), "1.5",
                      strand, ".", attributes]) + "\n"


CHAIN = [
    "#format: gtf\n",
    gtf_line("gene", 100, 900, 'gene_id "gene1"; gene_name "faster_gene1";'),
    gtf_line("transcript", 100, 900, 'gene_id "gene1"; transcript_id "tx1";'),
    gtf_line("exon", 100, 200, 'gene_id "gene1"; transcript_id "tx1"; exon_number "1";'),
    gtf_line("exon", 400, 500, 'gene_id "gene1"; transcript_id "tx1"; exon_number "2";'),
    gtf_line("exon", 800, 900, 'gene_id "gene1"; transcript_id "tx1"; exon_number "3";'),
    gtf_line("gene", 2000, 2100, 'gene_id "gene2";', strand="."),
    gtf_line("transcript", 2000, 2100, 'gene_id "gene2"; transcript_id "tx2";', strand="."),
    gtf_line("exon", 2000, 2100, 'gene_id "gene2"; transcript_id "tx2"; exon_number "1";',
             strand="."),
]


def records(lines):
    return [line.rstrip("\n").split("\t") for line in lines if not line.startswith("#")]


def test_a_three_exon_chain_becomes_three_single_exon_records():
    out = records(split_records(CHAIN))
    transcripts = [r for r in out if r[2] == "transcript"]
    exons = [r for r in out if r[2] == "exon"]
    assert len(transcripts) == len(exons) == 4
    assert [(r[3], r[4]) for r in transcripts] == [(r[3], r[4]) for r in exons]


def test_exon_coordinates_score_and_strand_are_unchanged():
    exons_in = [r for r in records(CHAIN) if r[2] == "exon"]
    exons_out = [r for r in records(split_records(CHAIN)) if r[2] == "exon"]
    assert [r[:8] for r in exons_out] == [r[:8] for r in exons_in]


def test_each_record_has_its_own_ids():
    out = records(split_records(CHAIN))
    transcript_ids = [r[8] for r in out if r[2] == "transcript"]
    assert len(set(transcript_ids)) == 4
    assert 'transcript_id "tx1.3"' in transcript_ids[2]
    assert 'gene_id "gene2.1"' in transcript_ids[3]


def test_header_lines_pass_through():
    assert list(split_records(CHAIN))[0] == "#format: gtf\n"
