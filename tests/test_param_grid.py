"""Unit tests for param_grid.py.

The grid decides which fastder runs exist and what each is called. An error
here changes the benchmark without failing, so the collapsing behaviour is
pinned below.
"""
import pytest

import param_grid as pg


class TestParamId:
    def test_values_use_the_abbreviation(self):
        assert pg.param_id({"min_coverage": 0.05, "min_length": 10}) == "mc0.05_ml10"

    def test_empty_combination_is_default(self):
        assert pg.param_id({}) == "default"

    def test_switch_records_both_states(self):
        assert pg.param_id({"no_stitch": False}) == "ns0"
        assert pg.param_id({"no_stitch": True}) == "ns1"

    def test_switch_states_get_separate_ids(self):
        on = pg.param_id({"min_coverage": 0.05, "no_stitch": True})
        off = pg.param_id({"min_coverage": 0.05, "no_stitch": False})
        assert on != off


class TestCliArgs:
    def test_values_are_passed_with_their_flag(self):
        assert pg.fastder_cli_args({"min_coverage": 0.05}) == "--min-coverage 0.05"

    def test_switch_is_passed_bare_when_true(self):
        assert pg.fastder_cli_args({"no_stitch": True}) == "--no-stitch"

    def test_switch_is_omitted_when_false(self):
        assert pg.fastder_cli_args({"no_stitch": False}) == ""

    def test_false_switch_leaves_no_stray_whitespace(self):
        args = pg.fastder_cli_args({"min_coverage": 0.05, "no_stitch": False})
        assert args == "--min-coverage 0.05"

    def test_min_junction_reads_is_a_value(self):
        assert pg.fastder_cli_args({"min_junction_reads": 5}) == "--min-junction-reads 5"


class TestDropIgnored:
    def test_no_stitch_drops_both_tolerances(self):
        combo = {"min_coverage": 0.05, "position_tolerance": 5,
                 "coverage_tolerance": 1000, "no_stitch": True}
        assert pg.drop_ignored(combo) == {"min_coverage": 0.05, "no_stitch": True}

    def test_tolerances_survive_when_stitching_is_on(self):
        combo = {"position_tolerance": 5, "no_stitch": False}
        assert pg.drop_ignored(combo) == combo

    def test_untouched_when_the_switch_is_absent(self):
        combo = {"min_coverage": 0.05, "position_tolerance": 5}
        assert pg.drop_ignored(combo) == combo


class TestBuildCombos:
    def test_empty_config_yields_one_default(self):
        assert pg.build_combos({}) == [{}]

    def test_cross_product_of_two_axes(self):
        combos = pg.build_combos({"min_coverage": [0.01, 0.05], "min_length": [10, 25]})
        assert len(combos) == 4
        assert {pg.param_id(c) for c in combos} == {
            "mc0.01_ml10", "mc0.01_ml25", "mc0.05_ml10", "mc0.05_ml25"}

    def test_unknown_keys_are_ignored(self):
        combos = pg.build_combos({"min_coverage": [0.05], "chromosomes": ["chr21"]})
        assert combos == [{"min_coverage": 0.05}]

    def test_no_stitch_collapses_the_tolerance_axis(self):
        # Four position_tolerance values against two no_stitch states would be
        # eight combinations, but --no-stitch ignores the tolerance, so the
        # four switched-on ones are one run.
        combos = pg.build_combos({"position_tolerance": [0, 5, 10, 20],
                                  "no_stitch": [False, True]})
        ids = [pg.param_id(c) for c in combos]
        assert len(ids) == 5
        assert ids.count("ns1") == 1
        assert sorted(i for i in ids if i.endswith("ns0")) == [
            "pt0_ns0", "pt10_ns0", "pt20_ns0", "pt5_ns0"]

    def test_collapsed_ids_omit_the_ignored_parameter(self):
        combos = pg.build_combos({"min_coverage": [0.05],
                                  "position_tolerance": [0, 5],
                                  "no_stitch": [False, True]})
        ids = {pg.param_id(c) for c in combos}
        assert "mc0.05_ns1" in ids
        assert not any(i.startswith("mc0.05_pt") and i.endswith("ns1") for i in ids)

    def test_every_combination_has_a_unique_id(self):
        combos = pg.build_combos({"min_coverage": [0.01, 0.05],
                                  "min_length": [10, 25],
                                  "position_tolerance": [0, 5, 10, 20],
                                  "no_stitch": [False, True]})
        ids = [pg.param_id(c) for c in combos]
        assert len(ids) == len(set(ids))
        assert len(ids) == 2 * 2 * 4 + 2 * 2
