"""The simulated reads are compressed, and the scenario step must keep them so.

open(path, "w") on a .gz name writes plain text, and STAR reads with zcat, so
the failure surfaces one rule later.
"""

import gzip

import make_scenario
from make_scenario import fastq_iter, filter_fastq, hardlink, passthrough

TEMPLATE = "ENST_TEMPLATE"
VARIANT = "ENST_VARIANT"


def record(transcript, index):
    return (f"@read{index}/{transcript};mate1:0-100;mate2:0-100\n"
            "ACGT\n"
            "+\n"
            "IIII\n")


def write_fastq(path, transcripts):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as fh:
        for i, transcript in enumerate(transcripts):
            fh.write(record(transcript, i))


def test_compressed_input_round_trips_compressed(tmp_path):
    src = tmp_path / "sample_01_1.fastq.gz"
    dst = tmp_path / "variant_only" / "sample_01_1.fastq.gz"
    dst.parent.mkdir()
    write_fastq(src, [TEMPLATE, VARIANT, TEMPLATE])

    filter_fastq(str(src), str(dst), {TEMPLATE})

    with gzip.open(dst, "rt") as fh:
        headers = [line for line in fh if line.startswith("@")]
    assert len(headers) == 1
    assert VARIANT in headers[0]


def test_plain_input_still_works(tmp_path):
    src = tmp_path / "sample_01_1.fastq"
    dst = tmp_path / "sample_01_1.filtered.fastq"
    write_fastq(src, [TEMPLATE, VARIANT])

    filter_fastq(str(src), str(dst), {TEMPLATE})

    assert dst.read_text().count("@read") == 1


def test_reader_takes_either_form(tmp_path):
    plain = tmp_path / "plain.fastq"
    packed = tmp_path / "packed.fastq.gz"
    write_fastq(plain, [VARIANT])
    write_fastq(packed, [VARIANT])
    assert list(fastq_iter(str(plain))) == list(fastq_iter(str(packed)))


def test_a_hard_linked_file_outlives_its_source(tmp_path):
    src = tmp_path / "sample_01_1.fastq.gz"
    dst = tmp_path / "template_and_variant" / "sample_01_1.fastq.gz"
    write_fastq(src, [TEMPLATE, VARIANT])

    hardlink(str(src), str(dst))
    src.unlink()

    with gzip.open(dst, "rt") as fh:
        assert fh.read().count("@read") == 2


def test_passthrough_links_rather_than_copies(tmp_path):
    src = tmp_path / "sample_01_1.fastq.gz"
    dst = tmp_path / "template_and_variant" / "sample_01_1.fastq.gz"
    write_fastq(src, [TEMPLATE, VARIANT])

    passthrough(str(src), str(dst))

    assert dst.is_symlink()
    with gzip.open(dst, "rt") as fh:
        assert fh.read().count("@read") == 2


GFF = (
    "chr21\tsim\tgene\t1\t100\t.\t+\t.\tgene_id=G1\n"
    f"chr21\tsim\ttranscript\t1\t100\t.\t+\t.\tgene_id=G1;transcript_id={TEMPLATE};template=TRUE\n"
    f"chr21\tsim\ttranscript\t1\t100\t.\t+\t.\tgene_id=G1;transcript_id={VARIANT};template=FALSE\n"
)


def run_main(monkeypatch, *arguments):
    monkeypatch.setattr("sys.argv", ["make_scenario.py", *map(str, arguments)])
    make_scenario.main()


def test_rebuilding_reads_leaves_the_truth_file_alone(tmp_path, monkeypatch):
    gff_in = tmp_path / "splicing_variants.gff3"
    gff_out = tmp_path / "variant_only" / "splicing_variants.gff3"
    fq_in = tmp_path / "sample_01_1.fastq.gz"
    fq_out = tmp_path / "variant_only" / "sample_01_1.fastq.gz"
    gff_in.write_text(GFF)
    write_fastq(fq_in, [TEMPLATE, VARIANT])

    run_main(monkeypatch, "--scenario", "variant_only",
             "--gff-in", gff_in, "--gff-out", gff_out)
    assert TEMPLATE not in gff_out.read_text()
    assert not fq_out.exists()
    truth_before = gff_out.stat().st_mtime_ns

    run_main(monkeypatch, "--scenario", "variant_only", "--gff-in", gff_in,
             "--fq1-in", fq_in, "--fq1-out", fq_out)
    with gzip.open(fq_out, "rt") as fh:
        assert fh.read().count("@read") == 1
    assert gff_out.stat().st_mtime_ns == truth_before
