"""Tests for savitzky_golay_smoother."""

import os
import tempfile

import numpy as np
import pytest

pytest.importorskip("pyopenms")


def _create_profile_mzml(path: str, n_spectra: int = 5) -> None:
    """Create a synthetic mzML with noisy profile spectra."""
    import pyopenms as oms

    exp = oms.MSExperiment()
    for i in range(n_spectra):
        spec = oms.MSSpectrum()
        spec.setMSLevel(1)
        spec.setRT(10.0 + i)

        # Dense m/z grid for profile data
        mz = np.linspace(400.0, 600.0, 2000)
        # Signal: two Gaussian peaks
        signal = (
            800.0 * np.exp(-0.5 * ((mz - 480.0) / 0.3) ** 2)
            + 1200.0 * np.exp(-0.5 * ((mz - 520.0) / 0.4) ** 2)
        )
        rng = np.random.default_rng(seed=100 + i)
        noise = rng.normal(0, 40.0, size=len(mz))
        intensity = np.maximum(signal + noise, 0.0).astype(np.float32)

        spec.set_peaks((mz.astype(np.float64), intensity))
        exp.addSpectrum(spec)

    oms.MzMLFile().store(path, exp)


class TestSavitzkyGolaySmoother:
    def test_spectra_count_matches(self):
        """Output spectrum count should match input."""
        from savitzky_golay_smoother import smooth_experiment

        with tempfile.TemporaryDirectory() as tmp:
            inp = os.path.join(tmp, "input.mzML")
            out = os.path.join(tmp, "output.mzML")
            _create_profile_mzml(inp, n_spectra=4)

            count = smooth_experiment(inp, out)
            assert count == 4

    def test_peaks_preserved(self):
        """Major peaks should still be present after smoothing."""
        import pyopenms as oms
        from savitzky_golay_smoother import smooth_experiment

        with tempfile.TemporaryDirectory() as tmp:
            inp = os.path.join(tmp, "input.mzML")
            out = os.path.join(tmp, "smoothed.mzML")
            _create_profile_mzml(inp, n_spectra=1)

            smooth_experiment(inp, out)

            smooth = oms.MSExperiment()
            oms.MzMLFile().load(out, smooth)
            mz, intensities = smooth[0].get_peaks()

            # Find index of peak near 520 (the stronger peak)
            idx_520 = np.argmin(np.abs(mz - 520.0))
            # The peak should still be prominent
            assert intensities[idx_520] > 500.0

    def test_noise_reduced(self):
        """Smoothing should reduce high-frequency noise."""
        import pyopenms as oms
        from savitzky_golay_smoother import smooth_experiment

        with tempfile.TemporaryDirectory() as tmp:
            inp = os.path.join(tmp, "input.mzML")
            out = os.path.join(tmp, "smoothed.mzML")
            _create_profile_mzml(inp)

            smooth_experiment(inp, out)

            orig = oms.MSExperiment()
            oms.MzMLFile().load(inp, orig)
            smooth = oms.MSExperiment()
            oms.MzMLFile().load(out, smooth)

            _, orig_int = orig[0].get_peaks()
            _, smooth_int = smooth[0].get_peaks()

            assert np.var(np.diff(smooth_int)) < np.var(np.diff(orig_int))

    def test_output_file_created(self):
        """Output mzML file should exist after smoothing."""
        from savitzky_golay_smoother import smooth_experiment

        with tempfile.TemporaryDirectory() as tmp:
            inp = os.path.join(tmp, "input.mzML")
            out = os.path.join(tmp, "output.mzML")
            _create_profile_mzml(inp, n_spectra=2)

            smooth_experiment(inp, out)
            assert os.path.isfile(out)

    def test_custom_parameters(self):
        """Custom frame_length and polynomial_order should be accepted."""
        from savitzky_golay_smoother import smooth_experiment

        with tempfile.TemporaryDirectory() as tmp:
            inp = os.path.join(tmp, "input.mzML")
            out = os.path.join(tmp, "output.mzML")
            _create_profile_mzml(inp, n_spectra=1)

            count = smooth_experiment(
                inp, out, frame_length=15, polynomial_order=4
            )
            assert count == 1
