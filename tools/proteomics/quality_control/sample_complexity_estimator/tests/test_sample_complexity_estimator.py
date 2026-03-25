"""Tests for sample_complexity_estimator."""

import json

import numpy as np
import pytest

pytest.importorskip("pyopenms")


class TestSampleComplexityEstimator:
    def _make_experiment(self, n_ms1=5, n_peaks=100):
        """Create a synthetic MSExperiment with MS1 spectra."""
        import pyopenms as oms

        exp = oms.MSExperiment()
        for i in range(n_ms1):
            spec = oms.MSSpectrum()
            spec.setMSLevel(1)
            spec.setRT(60.0 * i)
            mzs = np.array([100.0 + j * 1.5 for j in range(n_peaks)], dtype=np.float64)
            ints = np.array([1000.0 * (j + 1) for j in range(n_peaks)], dtype=np.float64)
            spec.set_peaks([mzs, ints])
            exp.addSpectrum(spec)
        return exp

    def test_estimate_complexity(self):
        from sample_complexity_estimator import estimate_complexity

        exp = self._make_experiment(n_ms1=3, n_peaks=50)
        result = estimate_complexity(exp)
        assert result["n_ms1_spectra"] == 3
        assert result["total_peaks"] == 150
        assert result["avg_peaks_per_spectrum"] == 50.0

    def test_complexity_score_low(self):
        from sample_complexity_estimator import estimate_complexity

        exp = self._make_experiment(n_ms1=2, n_peaks=80)
        result = estimate_complexity(exp)
        assert result["complexity_score"] == "low"

    def test_complexity_score_medium(self):
        from sample_complexity_estimator import estimate_complexity

        exp = self._make_experiment(n_ms1=2, n_peaks=500)
        result = estimate_complexity(exp)
        assert result["complexity_score"] == "medium"

    def test_intensity_threshold(self):
        from sample_complexity_estimator import estimate_complexity

        exp = self._make_experiment(n_ms1=1, n_peaks=10)
        # Threshold high enough to exclude some peaks
        result = estimate_complexity(exp, intensity_threshold=5000.0)
        assert result["total_peaks"] < 10

    def test_empty_experiment(self):
        import pyopenms as oms
        from sample_complexity_estimator import estimate_complexity

        exp = oms.MSExperiment()
        result = estimate_complexity(exp)
        assert result["n_ms1_spectra"] == 0
        assert result["complexity_score"] == "N/A"

    def test_per_spectrum_data(self):
        from sample_complexity_estimator import estimate_complexity

        exp = self._make_experiment(n_ms1=3, n_peaks=20)
        result = estimate_complexity(exp)
        assert len(result["per_spectrum"]) == 3
        assert all("rt" in s for s in result["per_spectrum"])
        assert all("n_peaks" in s for s in result["per_spectrum"])

    def test_write_json(self, tmp_path):
        from sample_complexity_estimator import estimate_complexity, write_json

        exp = self._make_experiment(n_ms1=2, n_peaks=30)
        result = estimate_complexity(exp)
        out = str(tmp_path / "complexity.json")
        write_json(result, out)
        with open(out) as fh:
            data = json.load(fh)
        assert data["n_ms1_spectra"] == 2
