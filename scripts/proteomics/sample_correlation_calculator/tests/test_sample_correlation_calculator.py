"""Tests for sample_correlation_calculator."""

import numpy as np
import pytest
from conftest import requires_pyopenms
from sample_correlation_calculator import compute_correlations, correlation_matrix


@requires_pyopenms
class TestSampleCorrelationCalculator:
    def _make_matrix(self):
        return np.array([
            [100.0, 200.0, 150.0],
            [300.0, 600.0, 450.0],
            [500.0, 1000.0, 750.0],
            [700.0, 1400.0, 1050.0],
        ])

    def test_pearson_self_correlation(self):
        matrix = self._make_matrix()
        col_names = ["s1", "s2", "s3"]
        results = compute_correlations(matrix, col_names, "pearson")
        self_corr = [r for r in results if r["sample_a"] == r["sample_b"]]
        for r in self_corr:
            assert abs(r["correlation"] - 1.0) < 1e-6

    def test_pearson_perfect_correlation(self):
        matrix = self._make_matrix()
        col_names = ["s1", "s2", "s3"]
        results = compute_correlations(matrix, col_names, "pearson")
        # s1 and s2 are perfectly linearly correlated (s2 = 2*s1)
        pair = next(r for r in results if r["sample_a"] == "s1" and r["sample_b"] == "s2")
        assert abs(pair["correlation"] - 1.0) < 1e-6

    def test_spearman(self):
        matrix = self._make_matrix()
        col_names = ["s1", "s2", "s3"]
        results = compute_correlations(matrix, col_names, "spearman")
        assert len(results) == 6  # 3 choose 2 + 3 diagonal

    def test_correlation_matrix_shape(self):
        matrix = self._make_matrix()
        col_names = ["s1", "s2", "s3"]
        corr_mat = correlation_matrix(matrix, col_names, "pearson")
        assert corr_mat.shape == (3, 3)
        # Diagonal should be 1
        np.testing.assert_allclose(np.diag(corr_mat), 1.0, atol=1e-6)

    def test_unknown_method(self):
        matrix = self._make_matrix()
        with pytest.raises(ValueError, match="Unknown method"):
            compute_correlations(matrix, ["s1", "s2", "s3"], "invalid")

    def test_with_nan(self):
        matrix = np.array([
            [100.0, 200.0],
            [np.nan, 400.0],
            [300.0, 600.0],
            [400.0, 800.0],
        ])
        results = compute_correlations(matrix, ["s1", "s2"], "pearson")
        # Should still compute using non-NaN rows
        pair = next(r for r in results if r["sample_a"] == "s1" and r["sample_b"] == "s2")
        assert not np.isnan(pair["correlation"])

    def test_pair_count(self):
        matrix = self._make_matrix()
        col_names = ["s1", "s2", "s3"]
        results = compute_correlations(matrix, col_names, "pearson")
        # n*(n+1)/2 = 6 pairs (including self)
        assert len(results) == 6
