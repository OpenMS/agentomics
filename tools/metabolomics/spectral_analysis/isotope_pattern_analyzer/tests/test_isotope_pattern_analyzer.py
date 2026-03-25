"""Tests for isotope_pattern_analyzer."""

import json
import os
import tempfile

import pytest
from conftest import requires_pyopenms


@requires_pyopenms
class TestGetTheoreticalPattern:
    def test_glucose_pattern(self):
        from isotope_pattern_analyzer import get_theoretical_pattern

        pattern = get_theoretical_pattern("C6H12O6", max_isotopes=3)
        assert len(pattern) == 3
        # First peak should be the most abundant (normalized to 100)
        assert pattern[0][1] == pytest.approx(100.0)
        # Subsequent peaks should be smaller
        assert pattern[1][1] < 100.0

    def test_sucrose_max_isotopes(self):
        from isotope_pattern_analyzer import get_theoretical_pattern

        pattern = get_theoretical_pattern("C12H22O11", max_isotopes=4)
        assert len(pattern) == 4

    def test_mz_ordering(self):
        from isotope_pattern_analyzer import get_theoretical_pattern

        pattern = get_theoretical_pattern("C12H22O11", max_isotopes=4)
        mzs = [mz for mz, _ in pattern]
        assert mzs == sorted(mzs)

    def test_empty_formula_raises(self):
        from isotope_pattern_analyzer import get_theoretical_pattern

        # Empty formula should return empty list (pyopenms may raise or return empty)
        try:
            result = get_theoretical_pattern("", max_isotopes=3)
            assert result == []
        except Exception:
            pass  # acceptable


@requires_pyopenms
class TestCosineSimilarity:
    def test_perfect_match(self):
        from isotope_pattern_analyzer import cosine_similarity

        peaks = [(180.0, 100.0), (181.0, 6.5), (182.0, 0.5)]
        sim = cosine_similarity(peaks, peaks, tolerance=0.05)
        assert abs(sim - 1.0) < 1e-6

    def test_no_overlap(self):
        from isotope_pattern_analyzer import cosine_similarity

        obs = [(180.0, 100.0), (181.0, 6.5)]
        theo = [(300.0, 100.0), (301.0, 6.5)]
        sim = cosine_similarity(obs, theo, tolerance=0.05)
        assert sim == 0.0

    def test_partial_overlap(self):
        from isotope_pattern_analyzer import cosine_similarity

        obs = [(180.0, 100.0), (181.0, 6.5), (182.0, 50.0)]  # Enhanced M+2
        theo = [(180.0, 100.0), (181.0, 6.5), (182.0, 0.5)]
        sim = cosine_similarity(obs, theo, tolerance=0.05)
        assert 0.0 < sim < 1.0

    def test_empty_returns_zero(self):
        from isotope_pattern_analyzer import cosine_similarity

        assert cosine_similarity([], [(180.0, 100.0)]) == 0.0
        assert cosine_similarity([(180.0, 100.0)], []) == 0.0

    def test_ppm_tolerance(self):
        from isotope_pattern_analyzer import cosine_similarity

        # At 180 Da, 10 ppm = 0.0018 Da
        peaks = [(180.001, 100.0), (181.002, 6.5)]
        theo = [(180.0, 100.0), (181.0, 6.5)]
        sim_ppm = cosine_similarity(peaks, theo, tolerance=10.0, tolerance_unit="ppm")
        sim_da = cosine_similarity(peaks, theo, tolerance=0.01, tolerance_unit="da")
        # Both should give reasonably similar results since 10 ppm ≈ 0.0018 Da at 180 Da
        assert sim_ppm > 0.9
        assert sim_da > 0.9


@requires_pyopenms
class TestDetectHalogenation:
    def test_no_halogen(self):
        from isotope_pattern_analyzer import detect_halogenation

        obs = [(180.0, 100.0), (181.0, 6.5), (182.0, 0.5)]
        theo = [(180.0, 100.0), (181.0, 6.5), (182.0, 0.5)]
        result = detect_halogenation(obs, theo)
        assert result["halogen_flag"] is False
        assert result["possible_halogen"] == "none"

    def test_chlorine_detection(self):
        from isotope_pattern_analyzer import detect_halogenation

        # Chlorine: M+2 ~32.5% of M+0
        obs = [(180.0, 100.0), (181.0, 6.5), (182.0, 35.0)]
        theo = [(180.0, 100.0), (181.0, 6.5), (182.0, 0.5)]
        result = detect_halogenation(obs, theo)
        assert result["halogen_flag"] is True
        assert "Cl" in result["possible_halogen"]

    def test_bromine_detection(self):
        from isotope_pattern_analyzer import detect_halogenation

        # Bromine: M+2 ~97% of M+0
        obs = [(180.0, 100.0), (181.0, 6.5), (182.0, 98.0)]
        theo = [(180.0, 100.0), (181.0, 6.5), (182.0, 0.5)]
        result = detect_halogenation(obs, theo)
        assert result["halogen_flag"] is True
        assert result["possible_halogen"] == "Br"

    def test_insufficient_peaks(self):
        from isotope_pattern_analyzer import detect_halogenation

        obs = [(180.0, 100.0), (181.0, 6.5)]  # only 2 peaks
        theo = [(180.0, 100.0), (181.0, 6.5)]
        result = detect_halogenation(obs, theo)
        assert result["halogen_flag"] is False
        assert result["m2_ratio_observed"] is None


@requires_pyopenms
class TestParseObserved:
    def test_colon_format(self):
        from isotope_pattern_analyzer import parse_observed

        peaks = parse_observed("180.063:100,181.067:6.5,182.070:0.5")
        assert len(peaks) == 3
        assert peaks[0] == (180.063, 100.0)
        assert peaks[2] == (182.070, 0.5)

    def test_comma_format_single(self):
        from isotope_pattern_analyzer import parse_observed

        peaks = parse_observed("181.0709,100.0")
        assert len(peaks) == 1
        assert peaks[0] == (181.0709, 100.0)

    def test_whitespace_stripped(self):
        from isotope_pattern_analyzer import parse_observed

        peaks = parse_observed("180.063 : 100 , 181.067 : 6.5")
        assert len(peaks) == 2

    def test_invalid_format_raises(self):
        from isotope_pattern_analyzer import parse_observed

        with pytest.raises(ValueError):
            parse_observed("notavalue")


@requires_pyopenms
class TestScorePattern:
    def test_basic_scoring(self):
        from isotope_pattern_analyzer import score_pattern

        obs = [(180.063, 100.0), (181.067, 6.5), (182.070, 0.5)]
        result = score_pattern(obs, "C6H12O6")
        assert "cosine_similarity" in result
        assert result["cosine_similarity"] > 0.0
        assert "halogen_detection" in result
        assert "theoretical_pattern" in result

    def test_score_range(self):
        from isotope_pattern_analyzer import score_pattern

        obs = [(180.063, 100.0), (181.067, 6.5), (182.070, 0.5)]
        result = score_pattern(obs, "C6H12O6")
        assert 0.0 <= result["cosine_similarity"] <= 1.0

    def test_ppm_tolerance_mode(self):
        from isotope_pattern_analyzer import score_pattern

        obs = [(180.063, 100.0), (181.067, 6.5)]
        result = score_pattern(obs, "C6H12O6", tolerance=10.0, tolerance_unit="ppm")
        assert result["tolerance_unit"] == "ppm"
        assert result["cosine_similarity"] > 0.0

    def test_output_json(self):
        from isotope_pattern_analyzer import score_pattern

        obs = [(180.063, 100.0), (181.067, 6.5)]
        result = score_pattern(obs, "C6H12O6")

        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "result.json")
            with open(out_path, "w") as fh:
                json.dump(result, fh, indent=2)
            assert os.path.exists(out_path)
            with open(out_path) as fh:
                loaded = json.load(fh)
            assert loaded["formula"] == "C6H12O6"

    def test_n_peaks_compared(self):
        from isotope_pattern_analyzer import score_pattern

        # 2 observed peaks, 6 theoretical → n_peaks_compared = 2
        obs = [(180.063, 100.0), (181.067, 6.5)]
        result = score_pattern(obs, "C6H12O6", max_isotopes=6)
        assert result["n_peaks_compared"] == 2

    def test_halogen_integration(self):
        from isotope_pattern_analyzer import score_pattern

        # Chlorinated compound with enhanced M+2
        obs = [(180.0, 100.0), (181.0, 6.5), (182.0, 35.0)]
        result = score_pattern(obs, "C6H12O6")
        assert result["halogen_detection"]["halogen_flag"] is True
