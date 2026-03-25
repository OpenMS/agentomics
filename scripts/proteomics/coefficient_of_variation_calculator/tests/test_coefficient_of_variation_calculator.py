"""Tests for coefficient_of_variation_calculator."""

import numpy as np
from coefficient_of_variation_calculator import calculate_cv, read_groups
from conftest import requires_pyopenms


@requires_pyopenms
class TestCoefficientOfVariationCalculator:
    def _make_data(self):
        matrix = np.array([
            [100.0, 110.0, 105.0, 200.0, 220.0, 210.0],
            [500.0, 500.0, 500.0, 1000.0, 1000.0, 1000.0],  # zero CV
        ])
        row_ids = ["prot1", "prot2"]
        col_names = ["s1", "s2", "s3", "s4", "s5", "s6"]
        groups = {"s1": "A", "s2": "A", "s3": "A", "s4": "B", "s5": "B", "s6": "B"}
        return matrix, row_ids, col_names, groups

    def test_basic_cv(self):
        matrix, row_ids, col_names, groups = self._make_data()
        results = calculate_cv(matrix, row_ids, col_names, groups)
        assert len(results) == 4  # 2 features x 2 groups

    def test_zero_cv(self):
        matrix, row_ids, col_names, groups = self._make_data()
        results = calculate_cv(matrix, row_ids, col_names, groups)
        prot2_a = next(r for r in results if r["feature"] == "prot2" and r["group"] == "A")
        assert abs(prot2_a["cv_percent"]) < 1e-6

    def test_cv_positive(self):
        matrix, row_ids, col_names, groups = self._make_data()
        results = calculate_cv(matrix, row_ids, col_names, groups)
        prot1_a = next(r for r in results if r["feature"] == "prot1" and r["group"] == "A")
        assert prot1_a["cv_percent"] > 0

    def test_n_values(self):
        matrix, row_ids, col_names, groups = self._make_data()
        results = calculate_cv(matrix, row_ids, col_names, groups)
        for r in results:
            assert r["n_values"] == 3

    def test_with_nan(self):
        matrix = np.array([[100.0, np.nan, 105.0, 200.0, 220.0, 210.0]])
        row_ids = ["prot1"]
        col_names = ["s1", "s2", "s3", "s4", "s5", "s6"]
        groups = {"s1": "A", "s2": "A", "s3": "A", "s4": "B", "s5": "B", "s6": "B"}
        results = calculate_cv(matrix, row_ids, col_names, groups)
        prot1_a = next(r for r in results if r["group"] == "A")
        assert prot1_a["n_values"] == 2

    def test_read_groups(self, tmp_path):
        gfile = str(tmp_path / "groups.tsv")
        with open(gfile, "w") as fh:
            fh.write("sample\tgroup\n")
            fh.write("s1\tA\n")
            fh.write("s2\tB\n")
        groups = read_groups(gfile)
        assert groups == {"s1": "A", "s2": "B"}

    def test_single_value_per_group(self):
        matrix = np.array([[100.0, 200.0]])
        row_ids = ["prot1"]
        col_names = ["s1", "s2"]
        groups = {"s1": "A", "s2": "B"}
        results = calculate_cv(matrix, row_ids, col_names, groups)
        # With only 1 value per group, CV should be NaN
        for r in results:
            assert np.isnan(r["cv_percent"])
