"""The simulated reads are stored compressed, so the scenario step has to keep
them that way.

The failure mode this guards is quiet: `open(path, "w")` on a name ending in
.gz writes plain text, and STAR is invoked with --readFilesCommand zcat, so the
alignment fails one rule later with nothing pointing back here.
"""

import gzip

from make_scenario import fastq_iter, filter_fastq, passthrough

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


def test_passthrough_links_rather_than_copies(tmp_path):
    src = tmp_path / "sample_01_1.fastq.gz"
    dst = tmp_path / "template_and_variant" / "sample_01_1.fastq.gz"
    write_fastq(src, [TEMPLATE, VARIANT])

    passthrough(str(src), str(dst))

    assert dst.is_symlink()
    with gzip.open(dst, "rt") as fh:
        assert fh.read().count("@read") == 2
