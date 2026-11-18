"""Tests for emit_lean_mm_rr.py: the MM file has per-sample read counts."""

import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "emit_lean_mm_rr.py"


def sj_line(chrom, start, end, unique, multi):
    return f"{chrom}\t{start}\t{end}\t1\t1\t0\t{unique}\t{multi}\t30\n"


def run(tmp_path, min_reads=1):
    a = tmp_path / "a.SJ.out.tab"
    b = tmp_path / "b.SJ.out.tab"
    a.write_text(sj_line("chr21", 100, 200, 7, 2) + sj_line("chr21", 300, 400, 1, 0)
                 + sj_line("chr1", 100, 200, 50, 0))
    b.write_text(sj_line("chr21", 100, 200, 12, 0))
    prefix = tmp_path / "junctions.ALL"
    subprocess.run([sys.executable, str(SCRIPT),
                    "--sample", "a", "--sj", str(a), "--sample", "b", "--sj", str(b),
                    "--chromosomes", "chr21", "--min-reads", str(min_reads),
                    "--out-prefix", str(prefix)], check=True, capture_output=True)
    lines = Path(str(prefix) + ".MM").read_text().splitlines()
    return lines[2], sorted(tuple(int(x) for x in line.split("\t")) for line in lines[3:])


def test_mm_has_unique_plus_multi_read_counts(tmp_path):
    header, triples = run(tmp_path)
    assert header == "2\t2\t3"
    assert triples == [(1, 1, 9), (1, 2, 12), (2, 1, 1)]


def test_min_reads_drops_low_pairs(tmp_path):
    _, triples = run(tmp_path, min_reads=2)
    assert triples == [(1, 1, 9), (1, 2, 12)]
