"""Tests for targeted_feature_extractor."""

from conftest import requires_pyopenms


@requires_pyopenms
class TestTargetedFeatureExtractor:
    def _make_experiment(self, target_mz=181.0707):
        """Create experiment with a peak at target_mz."""
        import numpy as np
        import pyopenms as oms

        exp = oms.MSExperiment()
        for i in range(5):
            spec = oms.MSSpectrum()
            spec.setMSLevel(1)
            spec.setRT(60.0 * i)
            mzs = np.array([target_mz - 1.0, target_mz, target_mz + 1.0], dtype=np.float64)
            ints = np.array([100.0, 5000.0 * (1 + i), 100.0], dtype=np.float64)
            spec.set_peaks([mzs, ints])
            exp.addSpectrum(spec)
        return exp

    def test_extract_eic(self):
        from targeted_feature_extractor import extract_eic

        exp = self._make_experiment(target_mz=181.0707)
        eic = extract_eic(exp, 181.0707, ppm=10.0)
        assert len(eic) == 5
        assert all(intensity > 0 for _, intensity in eic)

    def test_integrate_peak(self):
        from targeted_feature_extractor import integrate_peak

        eic = [(0.0, 100.0), (1.0, 100.0), (2.0, 100.0)]
        area = integrate_peak(eic)
        assert area == 200.0  # trapezoid: 2 intervals, 100 each

    def test_extract_targets(self):
        import pyopenms as oms
        from targeted_feature_extractor import extract_targets

        ef = oms.EmpiricalFormula("C6H12O6")
        target_mz = ef.getMonoWeight() + 1.007276
        exp = self._make_experiment(target_mz=target_mz)

        targets = [{"name": "Glucose", "formula": "C6H12O6"}]
        results = extract_targets(exp, targets, ppm=10.0)
        assert len(results) == 1
        assert results[0]["name"] == "Glucose"
        assert results[0]["peak_area"] > 0
        assert results[0]["max_intensity"] > 0

    def test_no_matching_peak(self):
        from targeted_feature_extractor import extract_eic

        exp = self._make_experiment(target_mz=181.0707)
        eic = extract_eic(exp, 500.0, ppm=5.0)
        assert all(intensity == 0.0 for _, intensity in eic)
