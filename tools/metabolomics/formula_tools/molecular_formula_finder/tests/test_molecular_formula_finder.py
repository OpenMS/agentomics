"""Tests for molecular_formula_finder."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestMolecularFormulaFinder:
    def test_find_glucose(self):
        from molecular_formula_finder import find_formulas

        constraints = {"C": (0, 12), "H": (0, 30), "N": (0, 5), "O": (0, 10)}
        results = find_formulas(180.0634, ppm=5.0, constraints=constraints)
        assert len(results) > 0
        formulas = [r["formula"] for r in results]
        assert "C6H12O6" in formulas

    def test_no_match(self):
        from molecular_formula_finder import find_formulas

        constraints = {"C": (0, 2), "H": (0, 2)}
        results = find_formulas(5000.0, ppm=1.0, constraints=constraints)
        assert len(results) == 0

    def test_parse_constraints(self):
        from molecular_formula_finder import parse_element_constraints

        constraints = parse_element_constraints("C:0-12,H:0-30,N:0-5")
        assert constraints["C"] == (0, 12)
        assert constraints["H"] == (0, 30)
        assert constraints["N"] == (0, 5)

    def test_senior_rule(self):
        from molecular_formula_finder import check_senior_rule

        assert check_senior_rule({"C": 6, "H": 12, "O": 6}) is True
        assert check_senior_rule({"H": 1}) is False

    def test_hc_ratio(self):
        from molecular_formula_finder import check_hc_ratio

        assert check_hc_ratio({"C": 6, "H": 12, "O": 6}) is True
        assert check_hc_ratio({"C": 1, "H": 100}) is False

    def test_result_keys(self):
        from molecular_formula_finder import find_formulas

        constraints = {"C": (6, 6), "H": (12, 12), "O": (6, 6)}
        results = find_formulas(180.0634, ppm=5.0, constraints=constraints)
        assert len(results) == 1
        r = results[0]
        assert "formula" in r
        assert "mass" in r
        assert "error_ppm" in r
        assert "passes_senior" in r

    def test_write_tsv(self):
        from molecular_formula_finder import write_tsv

        results = [{"formula": "C6H12O6", "mass": 180.0634, "error_ppm": 0.0,
                     "passes_senior": True, "passes_hc_ratio": True, "passes_nitrogen_rule": True}]
        with tempfile.TemporaryDirectory() as tmpdir:
            out = os.path.join(tmpdir, "formulas.tsv")
            write_tsv(results, out)
            assert os.path.exists(out)
