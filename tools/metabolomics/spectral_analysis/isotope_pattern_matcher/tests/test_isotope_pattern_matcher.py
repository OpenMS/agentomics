"""Tests for isotope_pattern_matcher."""

import pytest

pytest.importorskip("pyopenms")


class TestIsotopePatternMatcher:
    def test_glucose_pattern(self):
        from isotope_pattern_matcher import get_isotope_distribution

        dist = get_isotope_distribution("C6H12O6", max_isotopes=3)
        assert len(dist) == 3
        assert dist[0][1] == pytest.approx(100.0)
        assert dist[1][1] < dist[0][1]

    def test_pattern_mz_ordering(self):
        from isotope_pattern_matcher import get_isotope_distribution

        dist = get_isotope_distribution("C12H22O11", max_isotopes=4)
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

    def test_parse_peaks_invalid(self):
        from isotope_pattern_matcher import parse_peaks

        with pytest.raises(ValueError):
            parse_peaks(["181.0709"])
