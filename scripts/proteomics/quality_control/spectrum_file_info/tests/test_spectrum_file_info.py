"""Tests for spectrum_file_info."""

import pytest
from conftest import requires_pyopenms


@requires_pyopenms
class TestSpectrumFileInfo:
    def _make_experiment(self, n_spectra=5, ms_level=1):
        """Create a synthetic MSExperiment for testing."""
        import numpy as np
        import pyopenms as oms

        exp = oms.MSExperiment()
        for i in range(n_spectra):
            spec = oms.MSSpectrum()
            spec.setMSLevel(ms_level)
            spec.setRT(60.0 * i)
            mzs = np.array([100.0 + j for j in range(10)], dtype=np.float64)
            intensities = np.array([1000.0 * (j + 1) for j in range(10)], dtype=np.float64)
            spec.set_peaks([mzs, intensities])
            exp.addSpectrum(spec)
        return exp

    def test_summarise_nonempty(self):
        from spectrum_file_info import summarise_experiment

        exp = self._make_experiment(n_spectra=3)
        summary = summarise_experiment(exp)
        assert summary["n_spectra"] == 3
        assert 1 in summary["ms_levels"]

    def test_summarise_empty(self):
        import pyopenms as oms
        from spectrum_file_info import summarise_experiment

        exp = oms.MSExperiment()
        summary = summarise_experiment(exp)
        assert summary["n_spectra"] == 0

    def test_rt_range(self):
        from spectrum_file_info import summarise_experiment

        exp = self._make_experiment(n_spectra=5)
        summary = summarise_experiment(exp)
        rt_min, rt_max = summary["rt_range_sec"]
        assert rt_min == 0.0
        assert rt_max == 240.0

    def test_mz_range(self):
        from spectrum_file_info import summarise_experiment

        exp = self._make_experiment(n_spectra=2)
        summary = summarise_experiment(exp)
        mz_min, mz_max = summary["mz_range"]
        assert mz_min == pytest.approx(100.0)
        assert mz_max == pytest.approx(109.0)
