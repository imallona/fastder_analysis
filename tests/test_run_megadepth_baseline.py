"""Unit tests for run_megadepth_baseline.py.

Covers the pure-Python helpers (find_regions, mean_coverage, common_chroms)
and a CLI smoke test that builds two synthetic BigWigs in a tmp dir, runs
the script as a subprocess, and checks that the emitted GTF is structurally
valid and re-parseable.
"""
import os
import os.path as op
import subprocess
import sys

import numpy as np
import pyBigWig
import pytest

import run_megadepth_baseline as rmb


def _write_bw(path, chrom, length, intervals):
    """intervals: list of (start, end, value); zeros are implicit elsewhere."""
    bw = pyBigWig.open(path, "w")
    bw.addHeader([(chrom, length)])
    for start, end, value in intervals:
        bw.addEntries([chrom], [start], ends=[end], values=[float(value)])
    bw.close()


class TestFindRegions:
    def test_no_signal_returns_empty(self):
        cov = np.zeros(100, dtype=np.float64)
        assert rmb.find_regions(cov, cutoff=0.5, min_length=1) == []

    def test_single_run(self):
        cov = np.zeros(100, dtype=np.float64)
        cov[10:30] = 2.0
        regions = rmb.find_regions(cov, cutoff=1.0, min_length=1)
        assert regions == [(10, 30, 2.0)]

    def test_min_length_filters_short_runs(self):
        cov = np.zeros(100, dtype=np.float64)
        cov[10:14] = 2.0
        cov[40:60] = 2.0
        regions = rmb.find_regions(cov, cutoff=1.0, min_length=10)
        assert len(regions) == 1
        assert regions[0][:2] == (40, 60)

    def test_run_at_array_start(self):
        cov = np.zeros(50, dtype=np.float64)
        cov[0:10] = 3.0
        regions = rmb.find_regions(cov, cutoff=1.0, min_length=1)
        assert regions == [(0, 10, 3.0)]

    def test_run_at_array_end(self):
        cov = np.zeros(50, dtype=np.float64)
        cov[40:50] = 3.0
        regions = rmb.find_regions(cov, cutoff=1.0, min_length=1)
        assert regions == [(40, 50, 3.0)]

    def test_cutoff_is_inclusive(self):
        cov = np.full(20, 1.0, dtype=np.float64)
        regions = rmb.find_regions(cov, cutoff=1.0, min_length=1)
        assert regions == [(0, 20, 1.0)]
        regions = rmb.find_regions(cov, cutoff=1.5, min_length=1)
        assert regions == []

    def test_mean_coverage_value_returned(self):
        cov = np.zeros(40, dtype=np.float64)
        cov[10:20] = 1.0
        cov[20:30] = 5.0
        regions = rmb.find_regions(cov, cutoff=0.5, min_length=1)
        assert len(regions) == 1
        start, end, mean_cov = regions[0]
        assert (start, end) == (10, 30)
        assert mean_cov == pytest.approx(3.0)


class TestLibrarySize:
    """library_size = sum of (end - start) * value over the requested chroms."""

    def test_single_interval(self, tmp_path):
        bw = str(tmp_path / "a.all.bw")
        _write_bw(bw, "chr1", 100, [(10, 20, 4.0)])
        assert rmb.library_size(bw, ["chr1"]) == pytest.approx(40.0)

    def test_multiple_intervals_summed(self, tmp_path):
        bw = str(tmp_path / "a.all.bw")
        _write_bw(bw, "chr1", 200, [(10, 20, 4.0), (50, 70, 1.5)])
        assert rmb.library_size(bw, ["chr1"]) == pytest.approx(40.0 + 30.0)

    def test_only_requested_chroms_count(self, tmp_path):
        bw = str(tmp_path / "a.all.bw")
        bw_x = pyBigWig.open(bw, "w")
        bw_x.addHeader([("chr1", 100), ("chr2", 100)])
        bw_x.addEntries(["chr1"], [10], ends=[20], values=[4.0])
        bw_x.addEntries(["chr2"], [10], ends=[30], values=[2.0])
        bw_x.close()
        assert rmb.library_size(bw, ["chr1"]) == pytest.approx(40.0)
        assert rmb.library_size(bw, ["chr2"]) == pytest.approx(40.0)
        assert rmb.library_size(bw, ["chr1", "chr2"]) == pytest.approx(80.0)

    def test_missing_chrom_silently_skipped(self, tmp_path):
        bw = str(tmp_path / "a.all.bw")
        _write_bw(bw, "chr1", 100, [(10, 20, 4.0)])
        assert rmb.library_size(bw, ["chrZ"]) == 0.0


class TestMeanCpmCoverage:
    """Per-base mean across samples after per-sample CPM normalization."""

    def test_two_bigwigs_each_normalised_independently(self, tmp_path):
        # Each sample's library_size = (interval_length * value), so per-base
        # CPM inside the interval = value / (interval_length * value / 1e6)
        # = 1e6 / interval_length. Both samples have interval_length=10, so
        # per-base CPM is 1e5 for each, and the mean across samples is also
        # 1e5 inside the interval, 0 outside.
        chrom, length = "chr1", 100
        bw_a = str(tmp_path / "a.all.bw")
        bw_b = str(tmp_path / "b.all.bw")
        _write_bw(bw_a, chrom, length, [(10, 20, 4.0)])
        _write_bw(bw_b, chrom, length, [(10, 20, 2.0)])
        cpm_a = (10 * 4.0) / 1e6
        cpm_b = (10 * 2.0) / 1e6
        cov = rmb.mean_cpm_coverage([[bw_a], [bw_b]], chrom, length,
                                    [cpm_a, cpm_b])
        assert cov[10] == pytest.approx(1e5)
        assert cov[15] == pytest.approx(1e5)
        assert cov[5] == pytest.approx(0.0)
        assert cov[25] == pytest.approx(0.0)

    def test_zero_lib_size_skips_sample(self, tmp_path):
        chrom, length = "chr1", 50
        bw_a = str(tmp_path / "a.all.bw")
        bw_b = str(tmp_path / "b.all.bw")
        _write_bw(bw_a, chrom, length, [(10, 20, 5.0)])
        _write_bw(bw_b, chrom, length, [])
        cov = rmb.mean_cpm_coverage([[bw_a], [bw_b]], chrom, length,
                                    [50.0 * 1e-6, 0.0])
        # Sample b is skipped, sample a contributes 5.0 / cpm_factor over 10
        # bases. Then we still divide by len(bw_paths) = 2.
        assert cov[15] == pytest.approx((5.0 / (50.0 * 1e-6)) / 2)

    def test_missing_positions_count_as_zero(self, tmp_path):
        chrom, length = "chr1", 50
        bw_a = str(tmp_path / "a.all.bw")
        _write_bw(bw_a, chrom, length, [(10, 20, 5.0)])
        cov = rmb.mean_cpm_coverage([[bw_a]], chrom, length, [50.0 * 1e-6])
        # All positions outside the interval are 0
        assert cov[0] == pytest.approx(0.0)
        assert cov[40] == pytest.approx(0.0)


class TestLoadBigwigSamples:
    def test_finds_all_bw_first(self, tmp_path):
        (tmp_path / "s1.all.bw").touch()
        (tmp_path / "s1.plus.bw").touch()
        samples = rmb.load_bigwig_samples(str(tmp_path))
        assert samples == [[str(tmp_path / "s1.all.bw")]]

    def test_strand_tracks_group_into_one_sample(self, tmp_path):
        # A plus and a minus track are two strands of the same sample, so
        # they must come back as one sample, not two.
        (tmp_path / "s1.plus.bw").touch()
        (tmp_path / "s1.minus.bw").touch()
        samples = rmb.load_bigwig_samples(str(tmp_path))
        assert len(samples) == 1
        assert [op.basename(p) for p in samples[0]] == ["s1.minus.bw", "s1.plus.bw"]

    def test_no_bw_files_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            rmb.load_bigwig_samples(str(tmp_path))


class TestCommonChroms:
    def test_intersection(self, tmp_path):
        bw_a = str(tmp_path / "a.all.bw")
        bw_b = str(tmp_path / "b.all.bw")
        bw_x = pyBigWig.open(bw_a, "w")
        bw_x.addHeader([("chr1", 100), ("chr2", 100)])
        bw_x.close()
        bw_y = pyBigWig.open(bw_b, "w")
        bw_y.addHeader([("chr2", 100), ("chr3", 100)])
        bw_y.close()
        assert rmb.common_chroms([bw_a, bw_b]) == ["chr2"]


class TestCli:
    """End-to-end CLI test: build inputs, run as subprocess, parse output."""

    def test_emits_valid_gtf(self, tmp_path):
        chrom, length = "chr1", 200
        bw_a = str(tmp_path / "s1.all.bw")
        bw_b = str(tmp_path / "s2.all.bw")
        # Two clearly-expressed regions on both samples, well above 0.05 cutoff
        _write_bw(bw_a, chrom, length, [(20, 60, 5.0), (120, 180, 3.0)])
        _write_bw(bw_b, chrom, length, [(20, 60, 5.0), (120, 180, 3.0)])

        out_gtf = str(tmp_path / "output.gtf")
        script = op.join(op.dirname(op.dirname(op.abspath(__file__))),
                         "workflow", "scripts", "run_megadepth_baseline.py")
        subprocess.run(
            [sys.executable, script,
             "--bigwig-dir", str(tmp_path),
             "--out-gtf", out_gtf,
             "--cutoff", "0.05",
             "--min-length", "10"],
            check=True, capture_output=True,
        )

        with open(out_gtf) as fh:
            lines = [l for l in fh if not l.startswith("#")]
        # Each region emits gene + transcript + exon = 3 lines
        assert len(lines) == 6
        # All on the right chromosome and source tag
        for line in lines:
            fields = line.rstrip("\n").split("\t")
            assert len(fields) == 9
            assert fields[0] == chrom
            assert fields[1] == "megadepth_baseline"
            assert fields[2] in {"gene", "transcript", "exon"}
            assert fields[6] == "."  # unstranded
            assert "gene_id" in fields[8]
        # Check the two regions were emitted at the right coordinates (1-based GTF)
        gene_lines = [l for l in lines if l.split("\t")[2] == "gene"]
        coords = sorted((int(l.split("\t")[3]), int(l.split("\t")[4]))
                        for l in gene_lines)
        assert coords == [(21, 60), (121, 180)]

    def test_high_cpm_cutoff_emits_no_regions(self, tmp_path):
        # With CPM normalization, even tiny raw values inflate when divided
        # by a small library_size. To assert "no regions" we have to set the
        # cutoff above any per-base CPM the input could produce.
        chrom, length = "chr1", 100
        bw_a = str(tmp_path / "s1.all.bw")
        _write_bw(bw_a, chrom, length, [(10, 50, 1.0)])  # lib_size = 40

        out_gtf = str(tmp_path / "output.gtf")
        script = op.join(op.dirname(op.dirname(op.abspath(__file__))),
                         "workflow", "scripts", "run_megadepth_baseline.py")
        # Per-base CPM inside the interval is 1.0 / (40/1e6) = 25000. A
        # cutoff well above that gives no regions.
        subprocess.run(
            [sys.executable, script,
             "--bigwig-dir", str(tmp_path),
             "--out-gtf", out_gtf,
             "--cutoff", "1e9", "--min-length", "5"],
            check=True, capture_output=True,
        )
        with open(out_gtf) as fh:
            data_lines = [l for l in fh if not l.startswith("#")]
        assert data_lines == []

    def test_cpm_threshold_keeps_only_large_region(self, tmp_path):
        """Verify CPM scaling: a tiny region falls below cutoff while a
        100x larger region passes, with the cutoff set between their CPMs.
        """
        chrom, length = "chr1", 200
        bw = str(tmp_path / "s1.all.bw")
        _write_bw(bw, chrom, length, [(10, 50, 1.0), (60, 100, 100.0)])
        # lib_size = 40 * 1.0 + 40 * 100.0 = 4040; cpm_factor = 4040e-6
        # tiny region CPM ~ 247.5; big region CPM ~ 24752.5

        out_gtf = str(tmp_path / "output.gtf")
        script = op.join(op.dirname(op.dirname(op.abspath(__file__))),
                         "workflow", "scripts", "run_megadepth_baseline.py")
        subprocess.run(
            [sys.executable, script,
             "--bigwig-dir", str(tmp_path),
             "--out-gtf", out_gtf,
             "--cutoff", "1000", "--min-length", "5"],
            check=True, capture_output=True,
        )
        with open(out_gtf) as fh:
            gene_lines = [l for l in fh
                          if not l.startswith("#") and l.split("\t")[2] == "gene"]
        coords = sorted((int(l.split("\t")[3]), int(l.split("\t")[4]))
                        for l in gene_lines)
        # Only the big region (1-based GTF coordinates 61..100) survives
        assert coords == [(61, 100)]
