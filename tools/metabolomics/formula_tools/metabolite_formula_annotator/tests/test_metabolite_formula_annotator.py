"""Tests for metabolite_formula_annotator."""

import pytest

pytest.importorskip("pyopenms")


class TestMetaboliteFormulaAnnotator:
    def test_enumerate_formulas_glucose(self):
        from metabolite_formula_annotator import enumerate_formulas

        # Glucose: C6H12O6, mass ~180.0634
        candidates = enumerate_formulas(180.0634, ppm=5.0)
        formulas = [c["formula"] for c in candidates]
        assert "C6H12O6" in formulas

    def test_tolerance_filtering(self):
        from metabolite_formula_annotator import enumerate_formulas

        tight = enumerate_formulas(180.0634, ppm=1.0)
        loose = enumerate_formulas(180.0634, ppm=10.0)
        assert len(tight) <= len(loose)

    def test_score_isotope_fit(self):
        # Perfect match should score near 1.0
        import pyopenms as oms
        from metabolite_formula_annotator import score_isotope_fit

        ef = oms.EmpiricalFormula("C6H12O6")
        gen = oms.CoarseIsotopePatternGenerator(3)
        iso = ef.getIsotopeDistribution(gen)
        container = iso.getContainer()
        theo = [peak.getIntensity() for peak in container]
        max_t = max(theo)
        ratios = [t / max_t * 100.0 for t in theo]

        score = score_isotope_fit("C6H12O6", ratios)
        assert score > 0.99

    def test_annotate_features(self):
        import pyopenms as oms
        from metabolite_formula_annotator import PROTON, annotate_features

        ef = oms.EmpiricalFormula("C6H12O6")
        mz = ef.getMonoWeight() + PROTON

        features = [{"mz": str(mz), "rt": "60.0", "intensity": "1000"}]
        result = annotate_features(features, ppm=5.0, max_candidates=3)
        assert len(result) == 1
        assert len(result[0]["candidates"]) > 0
