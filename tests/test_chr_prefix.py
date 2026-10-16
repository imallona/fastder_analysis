import io

from chr_prefix import prefix_chromosomes, prefix_line


def test_prefix_line_adds_chr_once():
    assert prefix_line("19\tensembl\texon\t1\t9\n") == "chr19\tensembl\texon\t1\t9\n"
    assert prefix_line("chr19\tensembl\texon\t1\t9\n") == "chr19\tensembl\texon\t1\t9\n"


def test_prefix_line_keeps_comments_and_blank_lines():
    assert prefix_line("#!genome-build GRCh38\n") == "#!genome-build GRCh38\n"
    assert prefix_line("\n") == "\n"


def test_prefix_chromosomes_on_file_larger_than_read_buffer(tmp_path):
    record = "21\tensembl\texon\t100\t200\t.\t+\t.\tgene_id \"g\";\n"
    n_records = 2 * io.DEFAULT_BUFFER_SIZE // len(record) + 1
    src = tmp_path / "annotation.gtf"
    dst = tmp_path / "annotation.chr.gtf"
    src.write_text("#header\n" + record * n_records)
    assert src.stat().st_size > io.DEFAULT_BUFFER_SIZE

    prefix_chromosomes(src, dst)

    lines = dst.read_text().splitlines(keepends=True)
    assert lines[0] == "#header\n"
    assert lines[1:] == ["chr" + record] * n_records
