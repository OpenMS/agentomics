"""Tests for missing_value_imputation."""


import numpy as np
import pytest
from conftest import requires_pyopenms
from missing_value_imputation import (
    impute,
    impute_knn,
    impute_mindet,
    impute_minprob,
    read_matrix,
    write_matrix,
)


@requires_pyopenms
class TestMissingValueImputation:
    def _make_matrix_with_missing(self):
        matrix = np.array([
            [100.0, 200.0, np.nan],
            [150.0, np.nan, 300.0],
            [120.0, 180.0, 280.0],
            [np.nan, 210.0, 310.0],
        ])
        return matrix

    def test_mindet_no_nans(self):
        matrix = self._make_matrix_with_missing()
        result = impute_mindet(matrix)
        assert not np.any(np.isnan(result))

    def test_mindet_values(self):
        matrix = self._make_matrix_with_missing()
        result = impute_mindet(matrix)
        # Column 0 min = 100, Column 1 min = 180, Column 2 min = 280
        assert result[3, 0] == 100.0
        assert result[1, 1] == 180.0
        assert result[0, 2] == 280.0

    def test_minprob_no_nans(self):
        matrix = self._make_matrix_with_missing()
        result = impute_minprob(matrix)
        assert not np.any(np.isnan(result))

    def test_minprob_values_lower(self):
        matrix = self._make_matrix_with_missing()
        result = impute_minprob(matrix)
        # Imputed values should generally be lower than the mean
        col0_mean = np.nanmean(matrix[:, 0])
        assert result[3, 0] < col0_mean

    def test_knn_no_nans(self):
        matrix = self._make_matrix_with_missing()
        result = impute_knn(matrix, k=2)
        assert not np.any(np.isnan(result))

    def test_knn_reasonable_values(self):
        matrix = self._make_matrix_with_missing()
        result = impute_knn(matrix, k=2)
        # Imputed values should be within the range of observed values
        assert 50 < result[3, 0] < 500
        assert 50 < result[1, 1] < 500

    def test_impute_dispatch(self):
        matrix = self._make_matrix_with_missing()
        for method in ["mindet", "minprob", "knn"]:
            result = impute(matrix, method=method)
            assert not np.any(np.isnan(result))

    def test_unknown_method(self):
        matrix = self._make_matrix_with_missing()
        with pytest.raises(ValueError, match="Unknown imputation method"):
            impute(matrix, method="invalid")

    def test_read_write_roundtrip(self, tmp_path):
        row_ids = ["prot1", "prot2"]
        col_names = ["sample1", "sample2"]
        matrix = np.array([[100.0, 200.0], [300.0, 400.0]])
        outfile = str(tmp_path / "test.tsv")
        write_matrix(outfile, row_ids, col_names, matrix)
        r_ids, c_names, r_matrix = read_matrix(outfile)
        assert r_ids == row_ids
        assert c_names == col_names
        np.testing.assert_allclose(r_matrix, matrix, atol=0.01)

    def test_read_with_missing(self, tmp_path):
        outfile = str(tmp_path / "missing.tsv")
        with open(outfile, "w") as fh:
            fh.write("\tsample1\tsample2\n")
            fh.write("prot1\t100.0\tNA\n")
            fh.write("prot2\t\t200.0\n")
        _, _, matrix = read_matrix(outfile)
        assert np.isnan(matrix[0, 1])
        assert np.isnan(matrix[1, 0])

    def test_no_missing_unchanged(self):
        matrix = np.array([[1.0, 2.0], [3.0, 4.0]])
        result = impute_mindet(matrix)
        np.testing.assert_array_equal(result, matrix)
