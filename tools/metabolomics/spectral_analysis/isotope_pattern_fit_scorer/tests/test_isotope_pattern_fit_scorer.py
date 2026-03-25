"""Tests for isotope_pattern_fit_scorer."""

import json
import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestIsotopePatternFitScorer:
    def test_get_theoretical_pattern(self):
        from isotope_pattern_fit_scorer import get_theoretical_pattern

        pattern = get_theoretical_pattern("C6H12O6", max_isotopes=4)
        assert len(pattern) == 4
        # First peak should be the most abundant (normalized to 100)
        assert pattern[0][1] == pytest.approx(100.0)
        # Subsequent peaks should be smaller
        assert pattern[1][1] < 100.0

    def test_parse_observed(self):
        from isotope_pattern_fit_scorer import parse_observed

        peaks = parse_observed("180.063:100,181.067:6.5,182.070:0.5")
        assert len(peaks) == 3
        assert peaks[0] == (180.063, 100.0)
        assert peaks[1] == (181.067, 6.5)

    def test_parse_observed_with_spaces(self):
        from isotope_pattern_fit_scorer import parse_observed

        peaks = parse_observed("180.063:100 , 181.067:6.5")
        assert len(peaks) == 2

    def test_cosine_similarity_perfect(self):
        from isotope_pattern_fit_scorer import cosine_similarity

        peaks = [(180.0, 100.0), (181.0, 6.5), (182.0, 0.5)]
        sim = cosine_similarity(peaks, peaks, mz_tolerance=0.1)
        assert abs(sim - 1.0) < 1e-6

    def test_cosine_similarity_no_overlap(self):
        from isotope_pattern_fit_scorer import cosine_similarity

        obs = [(180.0, 100.0), (181.0, 6.5)]
        theo = [(300.0, 100.0), (301.0, 6.5)]
        sim = cosine_similarity(obs, theo, mz_tolerance=0.1)
        assert sim == 0.0

    def test_cosine_similarity_partial(self):
        from isotope_pattern_fit_scorer import cosine_similarity

        obs = [(180.0, 100.0), (181.0, 6.5), (182.0, 50.0)]  # Enhanced M+2
        theo = [(180.0, 100.0), (181.0, 6.5), (182.0, 0.5)]
        sim = cosine_similarity(obs, theo, mz_tolerance=0.1)
        assert 0.0 < sim < 1.0

    def test_detect_halogenation_no_halogen(self):
        from isotope_pattern_fit_scorer import detect_halogenation

        obs = [(180.0, 100.0), (181.0, 6.5), (182.0, 0.5)]
        theo = [(180.0, 100.0), (181.0, 6.5), (182.0, 0.5)]
        result = detect_halogenation(obs, theo)
        assert result["halogen_flag"] is False

    def test_detect_halogenation_chlorine(self):
        from isotope_pattern_fit_scorer import detect_halogenation

        obs = [(180.0, 100.0), (181.0, 6.5), (182.0, 35.0)]  # Enhanced M+2 (~35%)
        theo = [(180.0, 100.0), (181.0, 6.5), (182.0, 0.5)]   # Expected ~0.5%
        result = detect_halogenation(obs, theo)
        assert result["halogen_flag"] is True
        assert "Cl" in result["possible_halogen"]

    def test_detect_halogenation_bromine(self):
        from isotope_pattern_fit_scorer import detect_halogenation

        obs = [(180.0, 100.0), (181.0, 6.5), (182.0, 98.0)]  # Enhanced M+2 (~98%)
        theo = [(180.0, 100.0), (181.0, 6.5), (182.0, 0.5)]
        result = detect_halogenation(obs, theo)
        assert result["halogen_flag"] is True
        assert result["possible_halogen"] == "Br"

    def test_score_pattern(self):
        from isotope_pattern_fit_scorer import score_pattern

        result = score_pattern(
            "180.063:100,181.067:6.5,182.070:0.5",
            "C6H12O6",
        )
        assert "cosine_similarity" in result
        assert result["cosine_similarity"] > 0.0
        assert "halogen_detection" in result
        assert "observed_peaks" in result
        assert "theoretical_peaks" in result

    def test_score_pattern_output_json(self):
        from isotope_pattern_fit_scorer import score_pattern

        result = score_pattern("180.063:100,181.067:6.5", "C6H12O6")

        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "fit.json")
            with open(out_path, "w") as fh:
                json.dump(result, fh, indent=2)
            assert os.path.exists(out_path)
            with open(out_path) as fh:
                loaded = json.load(fh)
            assert loaded["formula"] == "C6H12O6"
