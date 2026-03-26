"""Tests for elution_peak_detector."""

import os
import tempfile

import pytest

pyopenms = pytest.importorskip("pyopenms")


class TestElutionPeakDetector:
    """Tests for elution peak detection functionality."""

    def test_create_synthetic_lcms_mzml(self):
        from elution_peak_detector import create_synthetic_lcms_mzml

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "run.mzML")
            create_synthetic_lcms_mzml(path, n_scans=20)

            exp = pyopenms.MSExperiment()
            pyopenms.MzMLFile().load(path, exp)
            assert exp.getNrSpectra() == 20
            assert exp.getSpectrum(0).getMSLevel() == 1

    def test_detect_elution_peaks_finds_peak(self):
        """Verify that a Gaussian elution profile at m/z=500 is detected as a peak."""
        from elution_peak_detector import create_synthetic_lcms_mzml, detect_elution_peaks

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "run.mzML")
            output_path = os.path.join(tmp, "peaks.mzML")

            create_synthetic_lcms_mzml(input_path, n_scans=20)

            peaks = detect_elution_peaks(
                input_path,
                output_path,
                width_filtering="off",
                mass_error_ppm=10.0,
                noise_threshold=100.0,
            )

            assert peaks >= 1
            assert os.path.exists(output_path)

    def test_detect_elution_peaks_output_chromatograms(self):
        """Verify that output mzML contains chromatograms for detected peaks."""
        from elution_peak_detector import create_synthetic_lcms_mzml, detect_elution_peaks

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "run.mzML")
            output_path = os.path.join(tmp, "peaks.mzML")

            create_synthetic_lcms_mzml(input_path, n_scans=20)

            peaks = detect_elution_peaks(
                input_path,
                output_path,
                width_filtering="off",
                mass_error_ppm=10.0,
                noise_threshold=100.0,
            )

            out_exp = pyopenms.MSExperiment()
            pyopenms.MzMLFile().load(output_path, out_exp)
            assert out_exp.getNrChromatograms() == peaks
            assert out_exp.getNrChromatograms() >= 1

    def test_detect_elution_peaks_centroid_mz(self):
        """Verify detected peak has centroid m/z near 500.0."""
        from elution_peak_detector import create_synthetic_lcms_mzml, detect_elution_peaks

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "run.mzML")
            output_path = os.path.join(tmp, "peaks.mzML")

            create_synthetic_lcms_mzml(input_path, n_scans=20)

            detect_elution_peaks(
                input_path,
                output_path,
                width_filtering="off",
                mass_error_ppm=10.0,
                noise_threshold=100.0,
            )

            out_exp = pyopenms.MSExperiment()
            pyopenms.MzMLFile().load(output_path, out_exp)
            chrom = out_exp.getChromatogram(0)
            product_mz = chrom.getProduct().getMZ()
            assert abs(product_mz - 500.0) < 1.0

    def test_detect_elution_peaks_high_noise_empty(self):
        """Verify that very high noise threshold returns zero peaks."""
        from elution_peak_detector import create_synthetic_lcms_mzml, detect_elution_peaks

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "run.mzML")
            output_path = os.path.join(tmp, "peaks.mzML")

            create_synthetic_lcms_mzml(input_path, n_scans=20)

            peaks = detect_elution_peaks(
                input_path,
                output_path,
                width_filtering="off",
                mass_error_ppm=10.0,
                noise_threshold=50000.0,
            )

            assert peaks == 0
