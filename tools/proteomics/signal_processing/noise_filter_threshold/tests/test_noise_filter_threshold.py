"""Tests for noise_filter_threshold."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


def _create_synthetic_mzml(path, mzs, intensities):
    """Create a synthetic mzML file with a single spectrum with given peaks."""
    import pyopenms as oms

    exp = oms.MSExperiment()
    spec = oms.MSSpectrum()
    spec.setMSLevel(1)
    spec.setRT(60.0)
    spec.set_peaks([mzs, intensities])

    exp.addSpectrum(spec)
    oms.MzMLFile().store(path, exp)


class TestFilterThreshold:
    def test_basic_threshold_filtering(self):
        """Peaks at intensities 50, 150, 500 with threshold=100; only 150 and 500 remain."""
        from noise_filter_threshold import filter_threshold

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.mzML")
            output_path = os.path.join(tmp, "output.mzML")

            mzs = [100.0, 200.0, 300.0]
            intensities = [50.0, 150.0, 500.0]
            _create_synthetic_mzml(input_path, mzs, intensities)

            count = filter_threshold(input_path, output_path, threshold=100.0)
            assert count == 1

            import pyopenms as oms

            exp = oms.MSExperiment()
            oms.MzMLFile().load(output_path, exp)

            spec = exp[0]
            _, filtered_intensities = spec.get_peaks()

            assert len(filtered_intensities) == 2
            for intensity in filtered_intensities:
                assert intensity >= 100.0

    def test_all_peaks_above_threshold(self):
        """When all peaks are above threshold, all remain."""
        from noise_filter_threshold import filter_threshold

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.mzML")
            output_path = os.path.join(tmp, "output.mzML")

            mzs = [100.0, 200.0, 300.0]
            intensities = [500.0, 600.0, 700.0]
            _create_synthetic_mzml(input_path, mzs, intensities)

            filter_threshold(input_path, output_path, threshold=100.0)

            import pyopenms as oms

            exp = oms.MSExperiment()
            oms.MzMLFile().load(output_path, exp)

            spec = exp[0]
            _, filtered_intensities = spec.get_peaks()
            assert len(filtered_intensities) == 3

    def test_all_peaks_below_threshold(self):
        """When all peaks are below threshold, none remain."""
        from noise_filter_threshold import filter_threshold

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.mzML")
            output_path = os.path.join(tmp, "output.mzML")

            mzs = [100.0, 200.0, 300.0]
            intensities = [10.0, 20.0, 30.0]
            _create_synthetic_mzml(input_path, mzs, intensities)

            filter_threshold(input_path, output_path, threshold=100.0)

            import pyopenms as oms

            exp = oms.MSExperiment()
            oms.MzMLFile().load(output_path, exp)

            spec = exp[0]
            _, filtered_intensities = spec.get_peaks()
            assert len(filtered_intensities) == 0

    def test_returns_spectrum_count(self):
        """Verify the function returns the number of spectra processed."""
        from noise_filter_threshold import filter_threshold

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.mzML")
            output_path = os.path.join(tmp, "output.mzML")

            mzs = [100.0, 200.0]
            intensities = [50.0, 150.0]
            _create_synthetic_mzml(input_path, mzs, intensities)

            count = filter_threshold(input_path, output_path, threshold=100.0)
            assert count == 1

    def test_output_file_created(self):
        """Verify the output file is created."""
        from noise_filter_threshold import filter_threshold

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.mzML")
            output_path = os.path.join(tmp, "output.mzML")

            mzs = [100.0]
            intensities = [50.0]
            _create_synthetic_mzml(input_path, mzs, intensities)

            filter_threshold(input_path, output_path, threshold=100.0)
            assert os.path.exists(output_path)
