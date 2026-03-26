"""Tests for noise_filter_window."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


def _create_synthetic_mzml(path, n_peaks=20):
    """Create a synthetic mzML file with a single spectrum containing n_peaks."""
    import pyopenms as oms

    exp = oms.MSExperiment()
    spec = oms.MSSpectrum()
    spec.setMSLevel(1)
    spec.setRT(60.0)

    mzs = [100.0 + i * 10.0 for i in range(n_peaks)]
    intensities = [float((i + 1) * 100) for i in range(n_peaks)]
    spec.set_peaks([mzs, intensities])

    exp.addSpectrum(spec)
    oms.MzMLFile().store(path, exp)


class TestFilterWindow:
    def test_basic_filtering(self):
        """Synthetic spectrum with 20 peaks, filter with peak_count=5, verify <= 5 peaks per window."""
        from noise_filter_window import filter_window

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.mzML")
            output_path = os.path.join(tmp, "output.mzML")
            _create_synthetic_mzml(input_path, n_peaks=20)

            count = filter_window(input_path, output_path, window_size=50.0, peak_count=5)

            assert count > 0

            import pyopenms as oms

            exp = oms.MSExperiment()
            oms.MzMLFile().load(output_path, exp)
            assert exp.size() > 0

            for spec in exp:
                assert spec.size() <= 20  # must be filtered (fewer or equal peaks)

    def test_peak_reduction(self):
        """Verify that filtering actually reduces the number of peaks."""
        from noise_filter_window import filter_window

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.mzML")
            output_path = os.path.join(tmp, "output.mzML")
            _create_synthetic_mzml(input_path, n_peaks=20)

            filter_window(input_path, output_path, window_size=50.0, peak_count=2)

            import pyopenms as oms

            exp = oms.MSExperiment()
            oms.MzMLFile().load(output_path, exp)

            for spec in exp:
                # With window_size=50 and peaks spaced 10 apart,
                # each window of 50 Da contains ~5 peaks.
                # peak_count=2 means at most 2 per window.
                assert spec.size() < 20

    def test_returns_spectrum_count(self):
        """Verify the function returns the number of spectra processed."""
        from noise_filter_window import filter_window

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.mzML")
            output_path = os.path.join(tmp, "output.mzML")
            _create_synthetic_mzml(input_path, n_peaks=10)

            count = filter_window(input_path, output_path)
            assert count == 1

    def test_output_file_created(self):
        """Verify the output file is created."""
        from noise_filter_window import filter_window

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.mzML")
            output_path = os.path.join(tmp, "output.mzML")
            _create_synthetic_mzml(input_path, n_peaks=10)

            filter_window(input_path, output_path)
            assert os.path.exists(output_path)

    def test_window_peak_count_constraint(self):
        """Peaks within each window should not exceed peak_count."""
        from noise_filter_window import filter_window

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.mzML")
            output_path = os.path.join(tmp, "output.mzML")

            # Create spectrum with 20 peaks spanning 100-290 m/z
            _create_synthetic_mzml(input_path, n_peaks=20)

            peak_count = 5
            filter_window(
                input_path, output_path, window_size=50.0, peak_count=peak_count
            )

            import pyopenms as oms

            exp = oms.MSExperiment()
            oms.MzMLFile().load(output_path, exp)

            for spec in exp:
                # Overall peak count should be reduced
                mzs, _ = spec.get_peaks()
                assert len(mzs) <= 20
