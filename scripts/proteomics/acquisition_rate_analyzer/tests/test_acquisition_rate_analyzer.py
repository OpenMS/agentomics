"""Tests for acquisition_rate_analyzer."""

from conftest import requires_pyopenms


@requires_pyopenms
class TestAcquisitionRateAnalyzer:
    def _make_experiment(self):
        import numpy as np
        import pyopenms as oms

        exp = oms.MSExperiment()
        rt = 0.0
        for i in range(5):
            ms1 = oms.MSSpectrum()
            ms1.setMSLevel(1)
            ms1.setRT(rt)
            mzs = np.array([100.0], dtype=np.float64)
            ints = np.array([1000.0], dtype=np.float64)
            ms1.set_peaks([mzs, ints])
            exp.addSpectrum(ms1)
            rt += 1.0

            for _ in range(3):
                ms2 = oms.MSSpectrum()
                ms2.setMSLevel(2)
                ms2.setRT(rt)
                ms2.set_peaks([np.array([200.0], dtype=np.float64),
                               np.array([500.0], dtype=np.float64)])
                exp.addSpectrum(ms2)
                rt += 0.3

        return exp

    def test_total_scans(self):
        from acquisition_rate_analyzer import analyze_acquisition_rates

        exp = self._make_experiment()
        result = analyze_acquisition_rates(exp)
        assert result["summary"]["total_scans"] == 20  # 5 MS1 + 15 MS2

    def test_ms1_ms2_counts(self):
        from acquisition_rate_analyzer import analyze_acquisition_rates

        exp = self._make_experiment()
        result = analyze_acquisition_rates(exp)
        assert result["summary"]["ms1_count"] == 5
        assert result["summary"]["ms2_count"] == 15

    def test_scan_records(self):
        from acquisition_rate_analyzer import analyze_acquisition_rates

        exp = self._make_experiment()
        result = analyze_acquisition_rates(exp)
        assert len(result["scans"]) == 20

    def test_empty_experiment(self):
        import pyopenms as oms
        from acquisition_rate_analyzer import analyze_acquisition_rates

        exp = oms.MSExperiment()
        result = analyze_acquisition_rates(exp)
        assert result["summary"]["total_scans"] == 0

    def test_ms2_per_cycle(self):
        from acquisition_rate_analyzer import analyze_acquisition_rates

        exp = self._make_experiment()
        result = analyze_acquisition_rates(exp)
        assert result["summary"]["avg_ms2_per_cycle"] == 3.0
