"""Tests for noise_filter_nlargest."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


def _create_synthetic_mzml(path, n_peaks=200):
    """Create a synthetic mzML file with a single spectrum containing n_peaks."""
    import pyopenms as oms

    exp = oms.MSExperiment()
    spec = oms.MSSpectrum()
    spec.setMSLevel(1)
    spec.setRT(60.0)

    mzs = [100.0 + i * 1.0 for i in range(n_peaks)]
    intensities = [float(i + 1) for i in range(n_peaks)]
    spec.set_peaks([mzs, intensities])

    exp.addSpectrum(spec)
    oms.MzMLFile().store(path, exp)


class TestFilterNLargest:
    def test_basic_nlargest_filtering(self):
        """Spectrum with 200 peaks, n=50, verify output has <= 50 peaks per spectrum."""
        from noise_filter_nlargest import filter_nlargest

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.mzML")
            output_path = os.path.join(tmp, "output.mzML")
            _create_synthetic_mzml(input_path, n_peaks=200)

            count = filter_nlargest(input_path, output_path, n=50)
            assert count == 1

            import pyopenms as oms

            exp = oms.MSExperiment()
            oms.MzMLFile().load(output_path, exp)

            for spec in exp:
                assert spec.size() <= 50

    def test_keeps_largest_peaks(self):
        """Verify that the largest peaks are kept."""
        from noise_filter_nlargest import filter_nlargest

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.mzML")
            output_path = os.path.join(tmp, "output.mzML")
            _create_synthetic_mzml(input_path, n_peaks=200)

            filter_nlargest(input_path, output_path, n=50)

            import pyopenms as oms

            exp = oms.MSExperiment()
            oms.MzMLFile().load(output_path, exp)

            spec = exp[0]
            _, intensities = spec.get_peaks()

            # All kept peaks should be among the 50 largest (intensity >= 151)
            for intensity in intensities:
                assert intensity >= 151.0

    def test_n_larger_than_peaks(self):
        """When n is larger than the number of peaks, all peaks are kept."""
        from noise_filter_nlargest import filter_nlargest

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.mzML")
            output_path = os.path.join(tmp, "output.mzML")
            _create_synthetic_mzml(input_path, n_peaks=10)

            filter_nlargest(input_path, output_path, n=100)

            import pyopenms as oms

            exp = oms.MSExperiment()
            oms.MzMLFile().load(output_path, exp)

            spec = exp[0]
            assert spec.size() == 10

    def test_returns_spectrum_count(self):
        """Verify the function returns the number of spectra processed."""
        from noise_filter_nlargest import filter_nlargest

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.mzML")
            output_path = os.path.join(tmp, "output.mzML")
            _create_synthetic_mzml(input_path, n_peaks=20)

            count = filter_nlargest(input_path, output_path, n=10)
            assert count == 1

    def test_output_file_created(self):
        """Verify the output file is created."""
        from noise_filter_nlargest import filter_nlargest

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.mzML")
            output_path = os.path.join(tmp, "output.mzML")
            _create_synthetic_mzml(input_path, n_peaks=20)

            filter_nlargest(input_path, output_path, n=10)
            assert os.path.exists(output_path)

    def test_multiple_spectra(self):
        """Verify filtering works across multiple spectra."""
        import pyopenms as oms
        from noise_filter_nlargest import filter_nlargest

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.mzML")
            output_path = os.path.join(tmp, "output.mzML")

            exp = oms.MSExperiment()
            for rt in [60.0, 120.0, 180.0]:
                spec = oms.MSSpectrum()
                spec.setMSLevel(1)
                spec.setRT(rt)
                mzs = [100.0 + i * 1.0 for i in range(100)]
                intensities = [float(i + 1) for i in range(100)]
                spec.set_peaks([mzs, intensities])
                exp.addSpectrum(spec)
            oms.MzMLFile().store(input_path, exp)

            count = filter_nlargest(input_path, output_path, n=10)
            assert count == 3

            exp_out = oms.MSExperiment()
            oms.MzMLFile().load(output_path, exp_out)
            for spec in exp_out:
                assert spec.size() <= 10
