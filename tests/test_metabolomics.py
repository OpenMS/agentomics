"""
Tests for metabolomics scripts.
Run with:  python -m pytest tests/ -v
"""

import sys
import os
import math

import pytest

# Allow importing from the scripts directory
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts", "metabolomics"))

try:
    import pyopenms as oms  # noqa: F401

    HAS_PYOPENMS = True
except ImportError:
    HAS_PYOPENMS = False

requires_pyopenms = pytest.mark.skipif(
    not HAS_PYOPENMS, reason="pyopenms not installed"
)


@requires_pyopenms
class TestMassAccuracyCalculator:
    def test_sequence_theoretical(self):
        from mass_accuracy_calculator import theoretical_mz_from_sequence

        mz = theoretical_mz_from_sequence("PEPTIDEK", 1)
        # The monoisotopic mass of PEPTIDEK is ~927.49 Da; +1 proton / 1
        assert 928.0 < mz < 929.0

    def test_formula_theoretical(self):
        from mass_accuracy_calculator import theoretical_mz_from_formula

        # Glucose C6H12O6 monoisotopic mass ~180.063 Da; +1 proton / 1 => ~181.070
        mz = theoretical_mz_from_formula("C6H12O6", 1)
        assert 181.0 < mz < 182.0

    def test_ppm_zero_error(self):
        from mass_accuracy_calculator import ppm_error

        theoretical = 500.0
        observed = 500.0
        assert ppm_error(theoretical, observed) == 0.0

    def test_ppm_positive_error(self):
        from mass_accuracy_calculator import ppm_error

        theoretical = 500.0
        observed = 500.001
        ppm = ppm_error(theoretical, observed)
        assert ppm > 0

    def test_ppm_negative_error(self):
        from mass_accuracy_calculator import ppm_error

        theoretical = 500.0
        observed = 499.999
        ppm = ppm_error(theoretical, observed)
        assert ppm < 0

    def test_ppm_known_value(self):
        from mass_accuracy_calculator import ppm_error

        # 1 ppm at 1000 Da => 0.001 Da shift
        theoretical = 1000.0
        observed = 1000.001
        ppm = ppm_error(theoretical, observed)
        assert abs(ppm - 1.0) < 0.001


@requires_pyopenms
class TestIsotopePatternMatcher:
    def test_glucose_pattern(self):
        from isotope_pattern_matcher import get_isotope_distribution

        dist = get_isotope_distribution("C6H12O6", max_isotopes=3)
        assert len(dist) == 3
        # M+0 should be the base peak
        assert dist[0][1] == pytest.approx(100.0)
        # Each successive peak must be less intense than the previous for small molecules
        assert dist[1][1] < dist[0][1]

    def test_pattern_mz_ordering(self):
        from isotope_pattern_matcher import get_isotope_distribution

        dist = get_isotope_distribution("C12H22O11", max_isotopes=4)  # sucrose
        mzs = [mz for mz, _ in dist]
        assert mzs == sorted(mzs)

    def test_cosine_similarity_perfect(self):
        from isotope_pattern_matcher import cosine_similarity

        peaks = [(100.0, 50.0), (101.0, 20.0), (102.0, 5.0)]
        sim = cosine_similarity(peaks, peaks, mz_tolerance=0.01)
        assert abs(sim - 1.0) < 1e-9

    def test_cosine_similarity_no_overlap(self):
        from isotope_pattern_matcher import cosine_similarity

        theoretical = [(100.0, 50.0), (101.0, 20.0)]
        observed = [(200.0, 50.0), (201.0, 20.0)]
        sim = cosine_similarity(theoretical, observed, mz_tolerance=0.01)
        assert sim == 0.0

    def test_parse_peaks(self):
        from isotope_pattern_matcher import parse_peaks

        result = parse_peaks(["181.0709,100.0", "182.0742,6.7"])
        assert len(result) == 2
        assert result[0] == (181.0709, 100.0)
        assert result[1] == (182.0742, 6.7)

    def test_parse_peaks_invalid(self):
        from isotope_pattern_matcher import parse_peaks

        with pytest.raises(ValueError):
            parse_peaks(["181.0709"])
