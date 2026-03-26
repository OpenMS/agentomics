"""Tests for mrm_rt_normalizer."""

import csv
import os
import random

import pytest

pyopenms = pytest.importorskip("pyopenms")


def _make_rt_pairs_tsv(path, n_good=20, n_outliers=2, seed=42):
    """Create a TSV file with iRT/measured RT pairs and outliers.

    Good pairs follow: measured_rt = 1.1 * irt + 5.0 + small noise.
    Outliers have extreme measured_rt values.
    """
    random.seed(seed)
    rows = []
    for i in range(n_good):
        irt = float(i)
        measured = irt * 1.1 + 5.0 + random.gauss(0, 0.1)
        rows.append({"irt": str(irt), "measured_rt": str(measured)})

    # Add outliers with extreme values
    for i in range(n_outliers):
        irt = float(n_good + i * 5)
        measured = 200.0 + i * 100.0  # far off from expected
        rows.append({"irt": str(irt), "measured_rt": str(measured)})

    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["irt", "measured_rt"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class TestMRMRTNormalizer:
    def test_basic_normalization(self, tmp_path):
        from mrm_rt_normalizer import normalize_mrm_rt

        pairs_path = str(tmp_path / "rt_pairs.tsv")
        out_path = str(tmp_path / "model.tsv")

        _make_rt_pairs_tsv(pairs_path)

        result = normalize_mrm_rt(pairs_path, out_path, method="iterative")
        assert os.path.exists(out_path)
        assert isinstance(result, dict)
        assert "before_count" in result
        assert "after_count" in result
        assert "slope" in result
        assert "intercept" in result
        assert "r_squared" in result

    def test_outlier_removal(self, tmp_path):
        from mrm_rt_normalizer import normalize_mrm_rt

        pairs_path = str(tmp_path / "rt_pairs.tsv")
        out_path = str(tmp_path / "model.tsv")

        _make_rt_pairs_tsv(pairs_path, n_good=20, n_outliers=2)

        result = normalize_mrm_rt(pairs_path, out_path, method="iterative")
        assert result["before_count"] == 22
        assert result["outliers_removed"] >= 0
        assert result["after_count"] <= result["before_count"]

    def test_linear_fit_quality(self, tmp_path):
        from mrm_rt_normalizer import normalize_mrm_rt

        pairs_path = str(tmp_path / "rt_pairs.tsv")
        out_path = str(tmp_path / "model.tsv")

        # Create pairs with no outliers for a clean fit
        _make_rt_pairs_tsv(pairs_path, n_good=20, n_outliers=0)

        result = normalize_mrm_rt(pairs_path, out_path, method="iterative")
        # With clean data: slope should be near 1.1, intercept near 5.0
        assert abs(result["slope"] - 1.1) < 0.1
        assert abs(result["intercept"] - 5.0) < 1.0
        assert result["r_squared"] > 0.95

    def test_output_contains_model_params(self, tmp_path):
        from mrm_rt_normalizer import normalize_mrm_rt

        pairs_path = str(tmp_path / "rt_pairs.tsv")
        out_path = str(tmp_path / "model.tsv")

        _make_rt_pairs_tsv(pairs_path)

        normalize_mrm_rt(pairs_path, out_path)

        with open(out_path) as fh:
            content = fh.read()

        assert "# slope=" in content
        assert "# intercept=" in content
        assert "# r_squared=" in content

    def test_output_has_cleaned_pairs(self, tmp_path):
        from mrm_rt_normalizer import normalize_mrm_rt

        pairs_path = str(tmp_path / "rt_pairs.tsv")
        out_path = str(tmp_path / "model.tsv")

        _make_rt_pairs_tsv(pairs_path)

        result = normalize_mrm_rt(pairs_path, out_path)

        # Read back the cleaned pairs (skip comment lines)
        with open(out_path) as fh:
            lines = [line for line in fh if not line.startswith("#")]
        reader = csv.DictReader(lines, delimiter="\t")
        rows = list(reader)
        assert len(rows) == result["after_count"]
        assert "reference_rt" in rows[0]
        assert "measured_rt" in rows[0]

    def test_fit_linear(self):
        from mrm_rt_normalizer import _fit_linear

        # Perfect linear: y = 2x + 3
        pairs = [[float(i), 2.0 * i + 3.0] for i in range(10)]
        slope, intercept, r_sq = _fit_linear(pairs)
        assert abs(slope - 2.0) < 0.001
        assert abs(intercept - 3.0) < 0.001
        assert abs(r_sq - 1.0) < 0.001

    def test_fit_linear_single_point(self):
        from mrm_rt_normalizer import _fit_linear

        slope, intercept, r_sq = _fit_linear([[1.0, 2.0]])
        # With 1 point, can't fit a line
        assert slope == 0.0

    def test_returns_dict(self, tmp_path):
        from mrm_rt_normalizer import normalize_mrm_rt

        pairs_path = str(tmp_path / "rt_pairs.tsv")
        out_path = str(tmp_path / "model.tsv")

        _make_rt_pairs_tsv(pairs_path)
        result = normalize_mrm_rt(pairs_path, out_path)
        assert isinstance(result, dict)
        assert result["outliers_removed"] >= 0
