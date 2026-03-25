"""Tests for quantification_normalizer."""

import numpy as np
import pytest
from conftest import requires_pyopenms
from quantification_normalizer import (
    normalize,
    normalize_median,
    normalize_quantile,
    normalize_total_intensity,
    read_matrix,
    write_matrix,
)


@requires_pyopenms
class TestQuantificationNormalizer:
    def _make_matrix(self):
        return np.array([
            [100.0, 200.0, 150.0],
            [300.0, 400.0, 350.0],
            [500.0, 600.0, 550.0],
            [700.0, 800.0, 750.0],
        ])

    def test_median_equal_medians(self):
        matrix = self._make_matrix()
        result = normalize_median(matrix)
        col_medians = np.median(result, axis=0)
        np.testing.assert_allclose(col_medians, col_medians[0], atol=1e-6)

    def test_quantile_equal_distributions(self):
        matrix = self._make_matrix()
        result = normalize_quantile(matrix)
        sorted_cols = np.sort(result, axis=0)
        for col in range(1, result.shape[1]):
            np.testing.assert_allclose(sorted_cols[:, 0], sorted_cols[:, col], atol=1e-6)

    def test_total_intensity_equal_sums(self):
        matrix = self._make_matrix()
        result = normalize_total_intensity(matrix)
        col_sums = np.sum(result, axis=0)
        np.testing.assert_allclose(col_sums, col_sums[0], atol=1e-6)

    def test_normalize_dispatch(self):
        matrix = self._make_matrix()
        for method in ["median", "quantile", "total_intensity"]:
            result = normalize(matrix, method=method)
            assert result.shape == matrix.shape

    def test_unknown_method(self):
        matrix = self._make_matrix()
        with pytest.raises(ValueError, match="Unknown normalization method"):
            normalize(matrix, method="invalid")

    def test_read_write_roundtrip(self, tmp_path):
        row_ids = ["prot1", "prot2"]
        col_names = ["s1", "s2"]
        matrix = np.array([[100.0, 200.0], [300.0, 400.0]])
        outfile = str(tmp_path / "test.tsv")
        write_matrix(outfile, row_ids, col_names, matrix)
        r_ids, c_names, r_matrix = read_matrix(outfile)
        assert r_ids == row_ids
        assert c_names == col_names
        np.testing.assert_allclose(r_matrix, matrix, atol=0.01)

    def test_preserves_shape(self):
        matrix = self._make_matrix()
        for method in ["median", "quantile", "total_intensity"]:
            assert normalize(matrix, method).shape == (4, 3)
