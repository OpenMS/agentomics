"""Tests for mass_decomposition_tool."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestMassDecompositionTool:
    def test_glucose_mass(self):
        from mass_decomposition_tool import decompose_mass

        # Glucose: C6H12O6 = 180.0634 Da
        results = decompose_mass(180.0634, tolerance=0.01)
        assert len(results) > 0
        formulas = [r["formula"] for r in results]
        assert "C6H12O6" in formulas

    def test_no_match(self):
        from mass_decomposition_tool import decompose_mass

        results = decompose_mass(5000.0, tolerance=0.001)
        assert len(results) == 0

    def test_custom_constraints(self):
        from mass_decomposition_tool import decompose_mass

        constraints = {"C": (6, 6), "H": (12, 12), "O": (6, 6)}
        results = decompose_mass(180.0634, tolerance=0.01, constraints=constraints)
        assert len(results) == 1
        assert results[0]["formula"] == "C6H12O6"

    def test_result_keys(self):
        from mass_decomposition_tool import decompose_mass

        results = decompose_mass(180.0634, tolerance=0.01)
        assert len(results) > 0
        for r in results:
            assert "formula" in r
            assert "mass" in r
            assert "error_da" in r

    def test_build_formula_string(self):
        from mass_decomposition_tool import _build_formula_string

        assert _build_formula_string({"C": 6, "H": 12, "O": 6}) == "C6H12O6"
        assert _build_formula_string({"C": 1, "H": 4}) == "CH4"
        assert _build_formula_string({}) == ""

    def test_write_tsv(self):
        from mass_decomposition_tool import write_tsv

        results = [{"formula": "C6H12O6", "mass": 180.0634, "error_da": 0.0}]
        with tempfile.TemporaryDirectory() as tmpdir:
            out = os.path.join(tmpdir, "decompositions.tsv")
            write_tsv(results, out)
            assert os.path.exists(out)
