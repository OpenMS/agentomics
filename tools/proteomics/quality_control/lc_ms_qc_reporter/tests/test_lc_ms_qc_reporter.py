"""Tests for lc_ms_qc_reporter."""

import pytest

pytest.importorskip("pyopenms")


class TestLcMsQcReporter:
    def _make_experiment(self):
        """Create a synthetic MSExperiment with MS1 and MS2 spectra."""
        import numpy as np
        import pyopenms as oms

        exp = oms.MSExperiment()
        for i in range(5):
            spec = oms.MSSpectrum()
            spec.setMSLevel(1)
            spec.setRT(60.0 * i)
            mzs = np.array([100.0 + j for j in range(10)], dtype=np.float64)
            intensities = np.array([1000.0 * (j + 1) for j in range(10)], dtype=np.float64)
            spec.set_peaks([mzs, intensities])
            exp.addSpectrum(spec)

            ms2 = oms.MSSpectrum()
            ms2.setMSLevel(2)
            ms2.setRT(60.0 * i + 1.0)
            prec = oms.Precursor()
            prec.setMZ(500.0 + i)
            prec.setCharge(2)
            ms2.setPrecursors([prec])
            mzs2 = np.array([200.0, 300.0, 400.0], dtype=np.float64)
            ints2 = np.array([500.0, 1000.0, 200.0], dtype=np.float64)
            ms2.set_peaks([mzs2, ints2])
            exp.addSpectrum(ms2)

        return exp

    def test_compute_qc_metrics(self):
        from lc_ms_qc_reporter import compute_qc_metrics

        exp = self._make_experiment()
        metrics = compute_qc_metrics(exp)
        assert metrics["ms1_count"] == 5
        assert metrics["ms2_count"] == 5
        assert metrics["total_spectra"] == 10

    def test_tic_stability(self):
        from lc_ms_qc_reporter import compute_qc_metrics

        exp = self._make_experiment()
        metrics = compute_qc_metrics(exp)
        # All MS1 spectra have same peaks, so CV should be 0
        assert metrics["tic_cv_percent"] == 0.0

    def test_charge_distribution(self):
        from lc_ms_qc_reporter import compute_qc_metrics

        exp = self._make_experiment()
        metrics = compute_qc_metrics(exp)
        assert "2" in metrics["charge_distribution"]
        assert metrics["charge_distribution"]["2"] == 5

    def test_rt_range(self):
        from lc_ms_qc_reporter import compute_qc_metrics

        exp = self._make_experiment()
        metrics = compute_qc_metrics(exp)
        assert metrics["rt_range_sec"][0] == 0.0
        assert metrics["rt_range_sec"][1] == 241.0

    def test_empty_experiment(self):
        import pyopenms as oms
        from lc_ms_qc_reporter import compute_qc_metrics

        exp = oms.MSExperiment()
        metrics = compute_qc_metrics(exp)
        assert metrics["ms1_count"] == 0
        assert metrics["ms2_count"] == 0
