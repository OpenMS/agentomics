"""Tests for run_comparison_reporter."""

import pytest

pytest.importorskip("pyopenms")


class TestRunComparisonReporter:
    def _make_experiment(self, rt_offset=0.0, prec_mz_base=500.0):
        import numpy as np
        import pyopenms as oms

        exp = oms.MSExperiment()
        for i in range(5):
            spec = oms.MSSpectrum()
            spec.setMSLevel(1)
            spec.setRT(60.0 * i + rt_offset)
            mzs = np.array([100.0 + j for j in range(10)], dtype=np.float64)
            ints = np.array([1000.0 * (j + 1) for j in range(10)], dtype=np.float64)
            spec.set_peaks([mzs, ints])
            exp.addSpectrum(spec)

            ms2 = oms.MSSpectrum()
            ms2.setMSLevel(2)
            ms2.setRT(60.0 * i + rt_offset + 1.0)
            prec = oms.Precursor()
            prec.setMZ(prec_mz_base + i)
            ms2.setPrecursors([prec])
            mzs2 = np.array([200.0], dtype=np.float64)
            ints2 = np.array([500.0], dtype=np.float64)
            ms2.set_peaks([mzs2, ints2])
            exp.addSpectrum(ms2)

        return exp

    def test_identical_runs(self):
        from run_comparison_reporter import compare_runs

        exp1 = self._make_experiment()
        exp2 = self._make_experiment()
        result = compare_runs(exp1, exp2)
        assert result["tic_correlation"] == 1.0
        assert result["shared_precursors"] == 5

    def test_different_precursors(self):
        from run_comparison_reporter import compare_runs

        exp1 = self._make_experiment(prec_mz_base=500.0)
        exp2 = self._make_experiment(prec_mz_base=600.0)
        result = compare_runs(exp1, exp2)
        assert result["shared_precursors"] == 0
        assert result["unique_to_run1"] == 5
        assert result["unique_to_run2"] == 5

    def test_rt_shift(self):
        from run_comparison_reporter import compare_runs

        exp1 = self._make_experiment(rt_offset=0.0)
        exp2 = self._make_experiment(rt_offset=30.0)
        result = compare_runs(exp1, exp2)
        assert result["rt_shift_sec"] == -30.0

    def test_pearson_correlation(self):
        from run_comparison_reporter import pearson_correlation

        assert abs(pearson_correlation([1, 2, 3], [1, 2, 3]) - 1.0) < 1e-6
        assert abs(pearson_correlation([1, 2, 3], [3, 2, 1]) + 1.0) < 1e-6
