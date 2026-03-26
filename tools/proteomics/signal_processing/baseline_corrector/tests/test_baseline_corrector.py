"""Tests for baseline_corrector."""

import os
import tempfile

import numpy as np
import pytest

pytest.importorskip("pyopenms")


def _create_baseline_mzml(path: str, n_spectra: int = 3) -> None:
    """Create synthetic mzML with spectra having a uniform baseline + peaks."""
    import pyopenms as oms

    exp = oms.MSExperiment()
    for i in range(n_spectra):
        spec = oms.MSSpectrum()
        spec.setMSLevel(1)
        spec.setRT(10.0 + i)

        # Dense m/z grid
        mz = np.linspace(400.0, 600.0, 2000)
        # Uniform baseline offset
        baseline = np.full_like(mz, 500.0)
        # Sharp peaks on top of baseline
        peak1 = 2000.0 * np.exp(-0.5 * ((mz - 470.0) / 0.2) ** 2)
        peak2 = 3000.0 * np.exp(-0.5 * ((mz - 530.0) / 0.3) ** 2)
        intensity = (baseline + peak1 + peak2).astype(np.float32)

        spec.set_peaks((mz.astype(np.float64), intensity))
        exp.addSpectrum(spec)

    oms.MzMLFile().store(path, exp)


class TestBaselineCorrector:
    def test_baseline_reduced(self):
        """Baseline intensity should be reduced after correction."""
        import pyopenms as oms
        from baseline_corrector import correct_baseline

        with tempfile.TemporaryDirectory() as tmp:
            inp = os.path.join(tmp, "input.mzML")
            out = os.path.join(tmp, "corrected.mzML")
            _create_baseline_mzml(inp)

            correct_baseline(inp, out, struct_element_length=3.0)

            corrected = oms.MSExperiment()
            oms.MzMLFile().load(out, corrected)
            mz, intensities = corrected[0].get_peaks()

            # Check a region far from peaks (around 450) -- baseline should
            # be much lower than the original 500
            idx_450 = np.argmin(np.abs(mz - 450.0))
            # Allow some tolerance: baseline should be reduced significantly
            assert intensities[idx_450] < 300.0

    def test_peaks_persist(self):
        """Peak regions should retain significant intensity after correction."""
        import pyopenms as oms
        from baseline_corrector import correct_baseline

        with tempfile.TemporaryDirectory() as tmp:
            inp = os.path.join(tmp, "input.mzML")
            out = os.path.join(tmp, "corrected.mzML")
            _create_baseline_mzml(inp)

            correct_baseline(inp, out)

            corrected = oms.MSExperiment()
            oms.MzMLFile().load(out, corrected)
            mz, intensities = corrected[0].get_peaks()

            # Peak at ~530 should still be prominent
            idx_530 = np.argmin(np.abs(mz - 530.0))
            assert intensities[idx_530] > 500.0

    def test_spectra_count(self):
        """Output should have the same number of spectra as input."""
        from baseline_corrector import correct_baseline

        with tempfile.TemporaryDirectory() as tmp:
            inp = os.path.join(tmp, "input.mzML")
            out = os.path.join(tmp, "corrected.mzML")
            _create_baseline_mzml(inp, n_spectra=5)

            count = correct_baseline(inp, out)
            assert count == 5

    def test_output_file_created(self):
        """Output mzML file should exist after correction."""
        from baseline_corrector import correct_baseline

        with tempfile.TemporaryDirectory() as tmp:
            inp = os.path.join(tmp, "input.mzML")
            out = os.path.join(tmp, "corrected.mzML")
            _create_baseline_mzml(inp, n_spectra=1)

            correct_baseline(inp, out)
            assert os.path.isfile(out)

    def test_custom_struct_element_length(self):
        """Custom structuring element length should be accepted."""
        from baseline_corrector import correct_baseline

        with tempfile.TemporaryDirectory() as tmp:
            inp = os.path.join(tmp, "input.mzML")
            out = os.path.join(tmp, "corrected.mzML")
            _create_baseline_mzml(inp, n_spectra=1)

            count = correct_baseline(inp, out, struct_element_length=5.0)
            assert count == 1
