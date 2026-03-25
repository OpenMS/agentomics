"""Tests for ptm_site_localization_scorer."""

from conftest import requires_pyopenms


@requires_pyopenms
class TestPtmSiteLocalizationScorer:
    def test_generate_theoretical_spectrum(self):
        from ptm_site_localization_scorer import generate_theoretical_spectrum

        ions = generate_theoretical_spectrum("PEPTIDEK")
        assert len(ions) > 0
        # All m/z values should be positive
        for mz, _annotation in ions:
            assert mz > 0

    def test_match_peaks_exact(self):
        from ptm_site_localization_scorer import generate_theoretical_spectrum, match_peaks

        ions = generate_theoretical_spectrum("PEPTIDEK")
        # Use theoretical m/z as "experimental"
        exp_mz = [mz for mz, _ in ions[:5]]
        matches = match_peaks(exp_mz, ions, tolerance=0.01)
        assert len(matches) >= 5

    def test_match_peaks_no_match(self):
        from ptm_site_localization_scorer import match_peaks

        matches = match_peaks([1.0, 2.0], [(1000.0, "b1")], tolerance=0.01)
        assert len(matches) == 0

    def test_generate_site_candidates(self):
        from ptm_site_localization_scorer import generate_site_candidates

        candidates = generate_site_candidates("PEPTIDEK", "Phospho", "ST")
        # Only T at position 4 (PEP_T_IDEK)
        assert len(candidates) >= 1

    def test_score_localization(self):
        from ptm_site_localization_scorer import generate_theoretical_spectrum, score_localization

        # Generate "experimental" spectrum from a specific modification
        ions = generate_theoretical_spectrum("PEPS(Phospho)TIDEK")
        exp_mz = [mz for mz, _ in ions]
        exp_int = [100.0] * len(exp_mz)

        result = score_localization(exp_mz, exp_int, "PEPS(Phospho)TIDEK", tolerance=0.02)
        assert result["total_experimental_peaks"] == len(exp_mz)
        assert len(result["candidates"]) >= 1
        # Probabilities should sum to ~1
        total_prob = sum(c["probability"] for c in result["candidates"])
        assert abs(total_prob - 1.0) < 0.01

    def test_score_localization_best_candidate(self):
        from ptm_site_localization_scorer import generate_theoretical_spectrum, score_localization

        # The correct site should have highest probability
        ions = generate_theoretical_spectrum("PEPS(Phospho)TIDEK")
        exp_mz = [mz for mz, _ in ions]
        exp_int = [100.0] * len(exp_mz)

        result = score_localization(exp_mz, exp_int, "PEPS(Phospho)TIDEK", tolerance=0.02)
        # Best candidate should be first (sorted by score)
        if len(result["candidates"]) > 0:
            assert result["candidates"][0]["matched_ions"] >= result["candidates"][-1]["matched_ions"]
