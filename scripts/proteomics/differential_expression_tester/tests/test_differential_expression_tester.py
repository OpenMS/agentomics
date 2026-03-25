"""Tests for differential_expression_tester."""

import math

import numpy as np
import pytest
from conftest import requires_pyopenms
from differential_expression_tester import (
    benjamini_hochberg,
    differential_expression,
    read_design,
)


@requires_pyopenms
class TestDifferentialExpressionTester:
    def _make_data(self):
        # 3 features, 6 samples (3 per condition)
        np.random.seed(42)
        matrix = np.array([
            [100, 110, 105, 200, 210, 205],   # clearly different
            [100, 105, 102, 103, 108, 101],    # not different
            [50, 55, 52, 500, 510, 505],       # very different
        ], dtype=float)
        row_ids = ["prot1", "prot2", "prot3"]
        col_names = ["s1", "s2", "s3", "s4", "s5", "s6"]
        design = {"s1": "A", "s2": "A", "s3": "A", "s4": "B", "s5": "B", "s6": "B"}
        return matrix, row_ids, col_names, design

    def test_basic_ttest(self):
        matrix, row_ids, col_names, design = self._make_data()
        results = differential_expression(matrix, row_ids, col_names, design, test="ttest")
        assert len(results) == 3
        assert all("pvalue" in r for r in results)
        assert all("log2fc" in r for r in results)
        assert all("adj_pvalue" in r for r in results)

    def test_significant_feature(self):
        matrix, row_ids, col_names, design = self._make_data()
        results = differential_expression(matrix, row_ids, col_names, design)
        # prot3 should be highly significant
        prot3 = next(r for r in results if r["feature"] == "prot3")
        assert prot3["pvalue"] < 0.01
        assert prot3["log2fc"] > 2.0  # ~10x increase

    def test_nonsignificant_feature(self):
        matrix, row_ids, col_names, design = self._make_data()
        results = differential_expression(matrix, row_ids, col_names, design)
        prot2 = next(r for r in results if r["feature"] == "prot2")
        assert prot2["pvalue"] > 0.05

    def test_welch_test(self):
        matrix, row_ids, col_names, design = self._make_data()
        results = differential_expression(matrix, row_ids, col_names, design, test="welch")
        assert len(results) == 3

    def test_bh_correction(self):
        pvalues = [0.01, 0.04, 0.03, 0.2]
        adjusted = benjamini_hochberg(pvalues)
        assert len(adjusted) == 4
        # Adjusted should be >= original
        for orig, adj in zip(pvalues, adjusted):
            assert adj >= orig or math.isnan(adj)

    def test_bh_all_nan(self):
        pvalues = [float("nan"), float("nan")]
        adjusted = benjamini_hochberg(pvalues)
        assert all(math.isnan(a) for a in adjusted)

    def test_wrong_conditions(self):
        matrix, row_ids, col_names, _ = self._make_data()
        design_3 = {"s1": "A", "s2": "B", "s3": "C", "s4": "A", "s5": "B", "s6": "C"}
        with pytest.raises(ValueError, match="Exactly 2 conditions"):
            differential_expression(matrix, row_ids, col_names, design_3)

    def test_design_file_read(self, tmp_path):
        design_file = str(tmp_path / "design.tsv")
        with open(design_file, "w") as fh:
            fh.write("sample\tcondition\n")
            fh.write("s1\tA\n")
            fh.write("s2\tB\n")
        design = read_design(design_file)
        assert design == {"s1": "A", "s2": "B"}
