"""Tests for peak_type_estimator."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestPeakTypeEstimator:
    def test_profile_spectrum_detected(self):
        import pyopenms as oms
        from peak_type_estimator import (
            create_synthetic_profile_spectrum,
            estimate_peak_type,
        )

        with tempfile.TemporaryDirectory() as tmp:
            mzml_path = os.path.join(tmp, "profile.mzML")
            out_path = os.path.join(tmp, "report.tsv")

            exp = oms.MSExperiment()
            exp.addSpectrum(create_synthetic_profile_spectrum())
            oms.MzMLFile().store(mzml_path, exp)

            result = estimate_peak_type(mzml_path, out_path)

            assert result["total_spectra"] == 1
            assert "profile" in result["type_counts"]
            assert result["type_counts"]["profile"] == 1

    def test_centroid_spectrum_detected(self):
        import pyopenms as oms
        from peak_type_estimator import (
            create_synthetic_centroid_spectrum,
            estimate_peak_type,
        )

        with tempfile.TemporaryDirectory() as tmp:
            mzml_path = os.path.join(tmp, "centroid.mzML")
            out_path = os.path.join(tmp, "report.tsv")

            exp = oms.MSExperiment()
            exp.addSpectrum(create_synthetic_centroid_spectrum())
            oms.MzMLFile().store(mzml_path, exp)

            result = estimate_peak_type(mzml_path, out_path)

            assert result["total_spectra"] == 1
            assert "centroid" in result["type_counts"]
            assert result["type_counts"]["centroid"] == 1

    def test_mixed_spectra(self):
        import pyopenms as oms
        from peak_type_estimator import (
            create_synthetic_centroid_spectrum,
            create_synthetic_profile_spectrum,
            estimate_peak_type,
        )

        with tempfile.TemporaryDirectory() as tmp:
            mzml_path = os.path.join(tmp, "mixed.mzML")
            out_path = os.path.join(tmp, "report.tsv")

            exp = oms.MSExperiment()
            exp.addSpectrum(create_synthetic_profile_spectrum())
            exp.addSpectrum(create_synthetic_centroid_spectrum())
            oms.MzMLFile().store(mzml_path, exp)

            result = estimate_peak_type(mzml_path, out_path)

            assert result["total_spectra"] == 2
            assert result["type_counts"].get("profile", 0) == 1
            assert result["type_counts"].get("centroid", 0) == 1

    def test_output_tsv_columns(self):
        import pyopenms as oms
        from peak_type_estimator import (
            create_synthetic_profile_spectrum,
            estimate_peak_type,
        )

        with tempfile.TemporaryDirectory() as tmp:
            mzml_path = os.path.join(tmp, "test.mzML")
            out_path = os.path.join(tmp, "report.tsv")

            exp = oms.MSExperiment()
            exp.addSpectrum(create_synthetic_profile_spectrum())
            oms.MzMLFile().store(mzml_path, exp)

            estimate_peak_type(mzml_path, out_path)

            with open(out_path) as f:
                header = f.readline().strip()
                assert header == "spectrum_index\tms_level\tpeak_type"
                lines = f.readlines()
                assert len(lines) == 1
