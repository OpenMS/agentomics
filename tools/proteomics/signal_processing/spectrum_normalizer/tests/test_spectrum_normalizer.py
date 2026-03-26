"""Tests for spectrum_normalizer."""

import os
import tempfile

import numpy as np
import pytest

pytest.importorskip("pyopenms")


def _create_test_mzml(path: str, n_spectra: int = 3) -> None:
    """Create a synthetic mzML with spectra having known intensities."""
    import pyopenms as oms

    exp = oms.MSExperiment()
    for i in range(n_spectra):
        spec = oms.MSSpectrum()
        spec.setMSLevel(1)
        spec.setRT(10.0 + i)

        mz = np.array([400.0, 450.0, 500.0, 550.0, 600.0])
        intensity = np.array(
            [100.0, 500.0, 1000.0, 300.0, 50.0], dtype=np.float32
        ) * (1.0 + 0.5 * i)

        spec.set_peaks((mz, intensity))
        exp.addSpectrum(spec)

    oms.MzMLFile().store(path, exp)


class TestSpectrumNormalizer:
    def test_to_one_max_intensity(self):
        """With method 'to_one', max intensity should be 1.0."""
        import pyopenms as oms
        from spectrum_normalizer import normalize_experiment

        with tempfile.TemporaryDirectory() as tmp:
            inp = os.path.join(tmp, "input.mzML")
            out = os.path.join(tmp, "normalized.mzML")
            _create_test_mzml(inp)

            normalize_experiment(inp, out, method="to_one")

            norm_exp = oms.MSExperiment()
            oms.MzMLFile().load(out, norm_exp)

            for i in range(norm_exp.size()):
                _, intensities = norm_exp[i].get_peaks()
                assert abs(max(intensities) - 1.0) < 1e-6

    def test_to_tic_sum_intensity(self):
        """With method 'to_TIC', sum of intensities should be ~1.0."""
        import pyopenms as oms
        from spectrum_normalizer import normalize_experiment

        with tempfile.TemporaryDirectory() as tmp:
            inp = os.path.join(tmp, "input.mzML")
            out = os.path.join(tmp, "normalized.mzML")
            _create_test_mzml(inp)

            normalize_experiment(inp, out, method="to_TIC")

            norm_exp = oms.MSExperiment()
            oms.MzMLFile().load(out, norm_exp)

            for i in range(norm_exp.size()):
                _, intensities = norm_exp[i].get_peaks()
                assert abs(sum(intensities) - 1.0) < 1e-4

    def test_spectra_count(self):
        """Output should have the same number of spectra."""
        from spectrum_normalizer import normalize_experiment

        with tempfile.TemporaryDirectory() as tmp:
            inp = os.path.join(tmp, "input.mzML")
            out = os.path.join(tmp, "normalized.mzML")
            _create_test_mzml(inp, n_spectra=4)

            count = normalize_experiment(inp, out)
            assert count == 4

    def test_output_file_created(self):
        """Output mzML file should exist after normalization."""
        from spectrum_normalizer import normalize_experiment

        with tempfile.TemporaryDirectory() as tmp:
            inp = os.path.join(tmp, "input.mzML")
            out = os.path.join(tmp, "normalized.mzML")
            _create_test_mzml(inp)

            normalize_experiment(inp, out)
            assert os.path.isfile(out)

    def test_relative_order_preserved(self):
        """Relative intensity ordering should be preserved after to_one."""
        import pyopenms as oms
        from spectrum_normalizer import normalize_experiment

        with tempfile.TemporaryDirectory() as tmp:
            inp = os.path.join(tmp, "input.mzML")
            out = os.path.join(tmp, "normalized.mzML")
            _create_test_mzml(inp, n_spectra=1)

            normalize_experiment(inp, out, method="to_one")

            norm_exp = oms.MSExperiment()
            oms.MzMLFile().load(out, norm_exp)
            _, intensities = norm_exp[0].get_peaks()

            # Original order: 100, 500, 1000, 300, 50
            # Peak at index 2 (1000) should be max, index 4 (50) should be min
            assert intensities[2] == max(intensities)
            assert intensities[4] == min(intensities)
