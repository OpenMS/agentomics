"""Tests for gaussian_smoother."""

import os
import tempfile

import numpy as np
import pytest

pytest.importorskip("pyopenms")


def _create_noisy_mzml(path: str, n_spectra: int = 5) -> None:
    """Create a synthetic mzML with noisy Gaussian profile spectra."""
    import pyopenms as oms

    exp = oms.MSExperiment()
    for i in range(n_spectra):
        spec = oms.MSSpectrum()
        spec.setMSLevel(1)
        spec.setRT(10.0 + i)

        # Dense m/z grid to simulate profile data
        mz = np.linspace(400.0, 600.0, 2000)
        # True signal: Gaussian peak centred at 500
        signal = 1000.0 * np.exp(-0.5 * ((mz - 500.0) / 0.5) ** 2)
        # Add random noise
        rng = np.random.default_rng(seed=42 + i)
        noise = rng.normal(0, 50.0, size=len(mz))
        intensity = np.maximum(signal + noise, 0.0).astype(np.float32)

        spec.set_peaks((mz.astype(np.float64), intensity))
        exp.addSpectrum(spec)

    oms.MzMLFile().store(path, exp)


class TestGaussianSmoother:
    def test_smooth_reduces_variance(self):
        """Smoothing should reduce intensity variance (noise)."""
        import pyopenms as oms
        from gaussian_smoother import smooth_experiment

        with tempfile.TemporaryDirectory() as tmp:
            inp = os.path.join(tmp, "noisy.mzML")
            out = os.path.join(tmp, "smoothed.mzML")
            _create_noisy_mzml(inp)

            count = smooth_experiment(inp, out, gaussian_width=0.2)
            assert count == 5

            # Load original and smoothed
            orig = oms.MSExperiment()
            oms.MzMLFile().load(inp, orig)
            smooth = oms.MSExperiment()
            oms.MzMLFile().load(out, smooth)

            # Compare variance of intensity differences between adjacent points
            for idx in range(orig.size()):
                _, orig_int = orig[idx].get_peaks()
                _, smooth_int = smooth[idx].get_peaks()
                orig_diff_var = np.var(np.diff(orig_int))
                smooth_diff_var = np.var(np.diff(smooth_int))
                assert smooth_diff_var < orig_diff_var

    def test_spectra_count_preserved(self):
        """Output should have the same number of spectra as input."""
        from gaussian_smoother import smooth_experiment

        with tempfile.TemporaryDirectory() as tmp:
            inp = os.path.join(tmp, "input.mzML")
            out = os.path.join(tmp, "output.mzML")
            _create_noisy_mzml(inp, n_spectra=3)

            count = smooth_experiment(inp, out)
            assert count == 3

    def test_output_file_created(self):
        """Output mzML file should exist after smoothing."""
        from gaussian_smoother import smooth_experiment

        with tempfile.TemporaryDirectory() as tmp:
            inp = os.path.join(tmp, "input.mzML")
            out = os.path.join(tmp, "output.mzML")
            _create_noisy_mzml(inp)

            smooth_experiment(inp, out)
            assert os.path.isfile(out)

    def test_custom_gaussian_width(self):
        """Wider Gaussian should produce smoother output."""
        import pyopenms as oms
        from gaussian_smoother import smooth_experiment

        with tempfile.TemporaryDirectory() as tmp:
            inp = os.path.join(tmp, "input.mzML")
            out_narrow = os.path.join(tmp, "narrow.mzML")
            out_wide = os.path.join(tmp, "wide.mzML")
            _create_noisy_mzml(inp)

            smooth_experiment(inp, out_narrow, gaussian_width=0.1)
            smooth_experiment(inp, out_wide, gaussian_width=1.0)

            narrow = oms.MSExperiment()
            oms.MzMLFile().load(out_narrow, narrow)
            wide = oms.MSExperiment()
            oms.MzMLFile().load(out_wide, wide)

            _, narrow_int = narrow[0].get_peaks()
            _, wide_int = wide[0].get_peaks()

            # Wider Gaussian should have lower point-to-point variance
            assert np.var(np.diff(wide_int)) <= np.var(np.diff(narrow_int))
