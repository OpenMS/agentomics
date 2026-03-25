"""Tests for feature_detection_proteomics."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestFeatureDetectionProteomics:
    def test_detect_features_returns_feature_map(self):
        import numpy as np
        import pyopenms as oms
        from feature_detection_proteomics import detect_features

        # Create a minimal synthetic experiment with a few peaks
        exp = oms.MSExperiment()
        for i in range(10):
            spec = oms.MSSpectrum()
            spec.setMSLevel(1)
            spec.setRT(60.0 + i * 2.0)
            mzs = np.array([500.0, 500.5, 501.0], dtype=np.float64)
            intensities = np.array([1e4, 5e3, 1e3], dtype=np.float64)
            spec.set_peaks([mzs, intensities])
            exp.addSpectrum(spec)

        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "test.mzML")
            output_path = os.path.join(tmpdir, "test.featureXML")
            oms.MzMLFile().store(input_path, exp)

            fm = detect_features(input_path, output_path)
            assert isinstance(fm, oms.FeatureMap)
            assert os.path.exists(output_path)
