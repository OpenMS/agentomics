"""Tests for metabolite_feature_detection."""

import os
import tempfile

from conftest import requires_pyopenms


@requires_pyopenms
class TestMetaboliteFeatureDetection:
    def test_detect_features_returns_feature_map(self):
        import numpy as np
        import pyopenms as oms
        from metabolite_feature_detection import detect_metabolite_features

        # Create a minimal synthetic experiment
        exp = oms.MSExperiment()
        for i in range(20):
            spec = oms.MSSpectrum()
            spec.setMSLevel(1)
            spec.setRT(30.0 + i * 3.0)
            mzs = np.array([180.063, 181.066, 182.070], dtype=np.float64)
            intensities = np.array([1e5, 1e4, 1e3], dtype=np.float64)
            spec.set_peaks([mzs, intensities])
            exp.addSpectrum(spec)

        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "test.mzML")
            output_path = os.path.join(tmpdir, "test.featureXML")
            oms.MzMLFile().store(input_path, exp)

            fm = detect_metabolite_features(input_path, output_path, noise_threshold=1e2)
            assert isinstance(fm, oms.FeatureMap)
            assert os.path.exists(output_path)
