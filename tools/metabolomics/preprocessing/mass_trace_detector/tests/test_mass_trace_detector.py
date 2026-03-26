"""Tests for mass_trace_detector."""

import os
import tempfile

import pytest

pyopenms = pytest.importorskip("pyopenms")


class TestMassTraceDetector:
    """Tests for mass trace detection functionality."""

    def test_create_synthetic_lcms_mzml(self):
        from mass_trace_detector import create_synthetic_lcms_mzml

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "run.mzML")
            create_synthetic_lcms_mzml(path, n_scans=20)

            exp = pyopenms.MSExperiment()
            pyopenms.MzMLFile().load(path, exp)
            assert exp.getNrSpectra() == 20
            assert exp.getSpectrum(0).getMSLevel() == 1

    def test_detect_mass_traces_finds_trace(self):
        """Verify that a persistent m/z=500 peak across 20 RT scans is detected as a trace."""
        from mass_trace_detector import create_synthetic_lcms_mzml, detect_mass_traces

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "run.mzML")
            output_path = os.path.join(tmp, "traces.mzML")

            create_synthetic_lcms_mzml(input_path, n_scans=20)

            traces = detect_mass_traces(
                input_path,
                output_path,
                mass_error_ppm=10.0,
                noise_threshold=100.0,
            )

            assert traces >= 1
            assert os.path.exists(output_path)

    def test_detect_mass_traces_output_has_chromatograms(self):
        """Verify that output mzML contains chromatograms for detected traces."""
        from mass_trace_detector import create_synthetic_lcms_mzml, detect_mass_traces

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "run.mzML")
            output_path = os.path.join(tmp, "traces.mzML")

            create_synthetic_lcms_mzml(input_path, n_scans=20)

            traces = detect_mass_traces(
                input_path,
                output_path,
                mass_error_ppm=10.0,
                noise_threshold=100.0,
            )

            # Load output and check chromatograms
            out_exp = pyopenms.MSExperiment()
            pyopenms.MzMLFile().load(output_path, out_exp)
            assert out_exp.getNrChromatograms() == traces
            assert out_exp.getNrChromatograms() >= 1

    def test_detect_mass_traces_centroid_mz(self):
        """Verify detected trace has centroid m/z near 500.0."""
        from mass_trace_detector import create_synthetic_lcms_mzml, detect_mass_traces

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "run.mzML")
            output_path = os.path.join(tmp, "traces.mzML")

            create_synthetic_lcms_mzml(input_path, n_scans=20)

            detect_mass_traces(
                input_path,
                output_path,
                mass_error_ppm=10.0,
                noise_threshold=100.0,
            )

            out_exp = pyopenms.MSExperiment()
            pyopenms.MzMLFile().load(output_path, out_exp)
            chrom = out_exp.getChromatogram(0)
            product_mz = chrom.getProduct().getMZ()
            assert abs(product_mz - 500.0) < 1.0

    def test_high_noise_threshold_filters_traces(self):
        """Verify that a very high noise threshold filters out weak traces."""
        from mass_trace_detector import create_synthetic_lcms_mzml, detect_mass_traces

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "run.mzML")
            output_path = os.path.join(tmp, "traces.mzML")

            create_synthetic_lcms_mzml(input_path, n_scans=20)

            # Use noise threshold higher than max intensity (10000)
            traces = detect_mass_traces(
                input_path,
                output_path,
                mass_error_ppm=10.0,
                noise_threshold=50000.0,
            )

            assert traces == 0
