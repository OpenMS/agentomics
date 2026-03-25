"""Tests for intensity_distribution_reporter."""

import numpy as np
from conftest import requires_pyopenms
from intensity_distribution_reporter import compute_intensity_stats, read_matrix


@requires_pyopenms
class TestIntensityDistributionReporter:
    def _make_matrix(self):
        return np.array([
            [100.0, 200.0, 150.0],
            [300.0, 400.0, 350.0],
            [500.0, 600.0, 550.0],
            [700.0, 800.0, 750.0],
        ])

    def test_basic_stats(self):
        matrix = self._make_matrix()
        col_names = ["s1", "s2", "s3"]
        stats = compute_intensity_stats(matrix, col_names)
        assert len(stats) == 3

    def test_mean_correct(self):
        matrix = self._make_matrix()
        col_names = ["s1", "s2", "s3"]
        stats = compute_intensity_stats(matrix, col_names)
        s1 = next(s for s in stats if s["sample"] == "s1")
        assert abs(s1["mean"] - 400.0) < 0.01  # mean of [100, 300, 500, 700]

    def test_min_max(self):
        matrix = self._make_matrix()
        col_names = ["s1", "s2", "s3"]
        stats = compute_intensity_stats(matrix, col_names)
        s1 = next(s for s in stats if s["sample"] == "s1")
        assert s1["min"] == 100.0
        assert s1["max"] == 700.0

    def test_n_values(self):
        matrix = self._make_matrix()
        col_names = ["s1", "s2", "s3"]
        stats = compute_intensity_stats(matrix, col_names)
        for s in stats:
            assert s["n_values"] == 4
            assert s["n_missing"] == 0

    def test_with_nan(self):
        matrix = np.array([
            [100.0, np.nan],
            [300.0, 400.0],
            [np.nan, 600.0],
        ])
        col_names = ["s1", "s2"]
        stats = compute_intensity_stats(matrix, col_names)
        s1 = next(s for s in stats if s["sample"] == "s1")
        assert s1["n_values"] == 2
        assert s1["n_missing"] == 1

    def test_all_nan_column(self):
        matrix = np.array([
            [np.nan, 200.0],
            [np.nan, 400.0],
        ])
        col_names = ["s1", "s2"]
        stats = compute_intensity_stats(matrix, col_names)
        s1 = next(s for s in stats if s["sample"] == "s1")
        assert s1["n_values"] == 0
        assert np.isnan(s1["mean"])

    def test_quartiles(self):
        matrix = np.array([[1.0], [2.0], [3.0], [4.0], [5.0], [6.0], [7.0], [8.0]])
        col_names = ["s1"]
        stats = compute_intensity_stats(matrix, col_names)
        s1 = stats[0]
        assert s1["q1"] == np.percentile([1, 2, 3, 4, 5, 6, 7, 8], 25)
        assert s1["q3"] == np.percentile([1, 2, 3, 4, 5, 6, 7, 8], 75)

    def test_read_matrix_roundtrip(self, tmp_path):
        outfile = str(tmp_path / "test.tsv")
        with open(outfile, "w") as fh:
            fh.write("\ts1\ts2\n")
            fh.write("p1\t100.0\t200.0\n")
            fh.write("p2\t300.0\tNA\n")
        row_ids, col_names, matrix = read_matrix(outfile)
        assert row_ids == ["p1", "p2"]
        assert col_names == ["s1", "s2"]
        assert np.isnan(matrix[1, 1])
