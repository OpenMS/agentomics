"""Tests for ms1_feature_intensity_tracker."""

import numpy as np
from conftest import requires_pyopenms


@requires_pyopenms
class TestMs1FeatureIntensityTracker:
    def _make_experiment(self, target_mz=500.0, target_rt=60.0, intensity=10000.0):
        """Create a synthetic MSExperiment with a known peak."""
        import pyopenms as oms

        exp = oms.MSExperiment()
        spec = oms.MSSpectrum()
        spec.setMSLevel(1)
        spec.setRT(target_rt)
        mzs = np.array([target_mz - 50, target_mz, target_mz + 50], dtype=np.float64)
        ints = np.array([100.0, intensity, 200.0], dtype=np.float64)
        spec.set_peaks([mzs, ints])
        exp.addSpectrum(spec)
        return exp

    def test_extract_intensity(self):
        from ms1_feature_intensity_tracker import extract_intensity

        exp = self._make_experiment(target_mz=500.0, intensity=5000.0)
        result = extract_intensity(exp, 500.0, ppm=10.0)
        assert result == 5000.0

    def test_extract_intensity_with_rt(self):
        from ms1_feature_intensity_tracker import extract_intensity

        exp = self._make_experiment(target_mz=500.0, target_rt=60.0, intensity=5000.0)
        # Within RT tolerance
        result = extract_intensity(exp, 500.0, ppm=10.0, target_rt=60.0, rt_tolerance=30.0)
        assert result == 5000.0
        # Outside RT tolerance
        result = extract_intensity(exp, 500.0, ppm=10.0, target_rt=200.0, rt_tolerance=30.0)
        assert result == 0.0

    def test_extract_intensity_not_found(self):
        from ms1_feature_intensity_tracker import extract_intensity

        exp = self._make_experiment(target_mz=500.0)
        result = extract_intensity(exp, 999.0, ppm=10.0)
        assert result == 0.0

    def test_track_features_from_experiments(self):
        from ms1_feature_intensity_tracker import track_features_from_experiments

        exp1 = self._make_experiment(target_mz=500.0, intensity=1000.0)
        exp2 = self._make_experiment(target_mz=500.0, intensity=2000.0)

        features = [{"feature_id": "F1", "mz": 500.0, "rt": None}]
        results = track_features_from_experiments(
            {"run1": exp1, "run2": exp2}, features, ppm=10.0,
        )
        assert len(results) == 1
        assert results[0]["run1"] == 1000.0
        assert results[0]["run2"] == 2000.0

    def test_load_features(self, tmp_path):
        from ms1_feature_intensity_tracker import load_features

        feat_path = str(tmp_path / "features.tsv")
        with open(feat_path, "w") as fh:
            fh.write("feature_id\tmz\trt\n")
            fh.write("F1\t500.0\t60.0\n")
            fh.write("F2\t600.0\t\n")

        features = load_features(feat_path)
        assert len(features) == 2
        assert features[0]["mz"] == 500.0
        assert features[0]["rt"] == 60.0
        assert features[1]["rt"] is None

    def test_write_tsv(self, tmp_path):
        from ms1_feature_intensity_tracker import write_tsv

        results = [{"feature_id": "F1", "mz": 500.0, "rt": 60.0, "run1": 1000.0, "run2": 2000.0}]
        out = str(tmp_path / "tracking.tsv")
        write_tsv(results, out)
        with open(out) as fh:
            lines = fh.readlines()
        assert len(lines) == 2
        assert "feature_id" in lines[0]

    def test_empty_experiment(self):
        import pyopenms as oms
        from ms1_feature_intensity_tracker import extract_intensity

        exp = oms.MSExperiment()
        result = extract_intensity(exp, 500.0, ppm=10.0)
        assert result == 0.0
