"""Unit tests for compute_library_sizes.py.

The library size is the whole-file sum of (end - start) * value, the same
quantity fastder accumulates and recount3 calls the AUC. Every tool in the
benchmark divides by these numbers, so a drift here shifts every threshold.
"""
import os.path as op
import subprocess
import sys

import pyBigWig
import pytest

import compute_library_sizes as cls


def _write_bw(path, chroms, entries):
    """chroms: list of (name, length). entries: list of (chrom, start, end, value)."""
    bw = pyBigWig.open(path, "w")
    bw.addHeader(list(chroms))
    for chrom, start, end, value in entries:
        bw.addEntries([chrom], [start], ends=[end], values=[float(value)])
    bw.close()


class TestLibrarySize:
    def test_single_interval(self, tmp_path):
        bw = str(tmp_path / "a.all.bw")
        _write_bw(bw, [("chr1", 100)], [("chr1", 10, 20, 4.0)])
        assert cls.library_size(bw) == pytest.approx(40.0)

    def test_base_weighted_not_a_plain_value_sum(self, tmp_path):
        # Two intervals of different lengths at the same value. A plain sum of
        # values would give 8.0; the base-weighted sum gives 120.0. bigWig's
        # sumData is documented only as "the sum of all values", so pin it.
        bw = str(tmp_path / "a.all.bw")
        _write_bw(bw, [("chr1", 200)], [("chr1", 0, 10, 4.0), ("chr1", 20, 40, 4.0)])
        assert cls.library_size(bw) == pytest.approx(40.0 + 80.0)

    def test_sums_across_chromosomes(self, tmp_path):
        bw = str(tmp_path / "a.all.bw")
        _write_bw(bw, [("chr1", 100), ("chr2", 100)],
                  [("chr1", 10, 20, 4.0), ("chr2", 10, 30, 2.0)])
        assert cls.library_size(bw) == pytest.approx(40.0 + 40.0)

    def test_empty_file_is_zero(self, tmp_path):
        bw = str(tmp_path / "a.all.bw")
        _write_bw(bw, [("chr1", 100)], [])
        assert cls.library_size(bw) == 0.0


class TestSampleId:
    def test_strips_strand_and_suffix(self):
        assert cls.sample_id("/x/SRR1.all.bw") == "SRR1"
        assert cls.sample_id("/x/SRR1.plus.bw") == "SRR1"
        assert cls.sample_id("/x/SRR1.minus.bw") == "SRR1"
        assert cls.sample_id("/x/SRR1.unique.bw") == "SRR1"

    def test_strips_bare_suffix(self):
        assert cls.sample_id("/x/SRR1.bw") == "SRR1"

    def test_stranded_tracks_share_a_sample_id(self):
        assert cls.sample_id("/x/S.plus.bw") == cls.sample_id("/x/S.minus.bw")


class TestCli:
    def test_writes_one_row_per_bigwig(self, tmp_path):
        _write_bw(str(tmp_path / "s1.all.bw"), [("chr1", 100)],
                  [("chr1", 10, 20, 4.0)])
        _write_bw(str(tmp_path / "s2.all.bw"), [("chr1", 100)],
                  [("chr1", 0, 50, 2.0)])
        out = tmp_path / "library_sizes.tsv"
        script = op.join(op.dirname(op.dirname(op.abspath(__file__))),
                         "workflow", "scripts", "compute_library_sizes.py")
        subprocess.run([sys.executable, script,
                        "--bigwig-dir", str(tmp_path),
                        "--out", str(out)],
                       check=True, capture_output=True)

        rows = [l.rstrip("\n").split("\t") for l in out.read_text().splitlines()]
        assert rows[0] == ["bigwig", "sample", "library_size"]
        sizes = {r[1]: float(r[2]) for r in rows[1:]}
        assert sizes == {"s1": pytest.approx(40.0), "s2": pytest.approx(100.0)}

    def test_no_bigwigs_exits_nonzero(self, tmp_path):
        script = op.join(op.dirname(op.dirname(op.abspath(__file__))),
                         "workflow", "scripts", "compute_library_sizes.py")
        result = subprocess.run([sys.executable, script,
                                 "--bigwig-dir", str(tmp_path),
                                 "--out", str(tmp_path / "out.tsv")],
                                capture_output=True)
        assert result.returncode != 0
