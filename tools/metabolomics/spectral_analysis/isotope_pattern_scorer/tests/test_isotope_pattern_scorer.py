"""Tests for isotope_pattern_scorer."""

from conftest import requires_pyopenms


@requires_pyopenms
class TestIsotopePatternScorer:
    def test_parse_observed(self):
        from isotope_pattern_scorer import parse_observed

        peaks = parse_observed("180.063:100,181.067:6.5,182.070:0.5")
        assert len(peaks) == 3
        assert peaks[0] == (180.063, 100.0)

    def test_get_theoretical_pattern(self):
        from isotope_pattern_scorer import get_theoretical_pattern

        pattern = get_theoretical_pattern("C6H12O6", n_peaks=3)
        assert len(pattern) == 3
        # First peak should be the most intense (normalized to 100)
        assert pattern[0][1] == 100.0

    def test_perfect_match_score(self):
        from isotope_pattern_scorer import get_theoretical_pattern, score_pattern

        theo = get_theoretical_pattern("C6H12O6", n_peaks=3)
        # Use theoretical as observed
        result = score_pattern(theo, theo)
        assert result["cosine_score"] > 0.999

    def test_score_range(self):
        from isotope_pattern_scorer import get_theoretical_pattern, score_pattern

        theo = get_theoretical_pattern("C6H12O6", n_peaks=3)
        observed = [(180.063, 100.0), (181.067, 50.0), (182.070, 30.0)]
        result = score_pattern(observed, theo)
        assert 0.0 <= result["cosine_score"] <= 1.0

    def test_empty_observed(self):
        from isotope_pattern_scorer import score_pattern

        result = score_pattern([], [])
        assert result["cosine_score"] == 0.0
        assert result["n_peaks_compared"] == 0

    def test_peaks_detail(self):
        from isotope_pattern_scorer import get_theoretical_pattern, score_pattern

        theo = get_theoretical_pattern("C6H12O6", n_peaks=2)
        observed = [(180.063, 100.0), (181.067, 7.0)]
        result = score_pattern(observed, theo)
        assert result["n_peaks_compared"] == 2
        assert len(result["peaks"]) == 2
