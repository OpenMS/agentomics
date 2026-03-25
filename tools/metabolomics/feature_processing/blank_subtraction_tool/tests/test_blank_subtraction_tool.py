"""Tests for blank_subtraction_tool."""

import pytest

pytest.importorskip("pyopenms")


class TestBlankSubtractionTool:
    def test_remove_blank_features(self):
        from blank_subtraction_tool import subtract_blanks

        sample = [
            {"mz": "100.0", "rt": "60.0", "intensity": "1000"},
            {"mz": "200.0", "rt": "120.0", "intensity": "500"},
        ]
        blank = [
            {"mz": "100.0", "rt": "60.0", "intensity": "800"},  # high in blank
        ]

        cleaned = subtract_blanks(sample, blank, fold_change=3.0)
        # Feature at 100.0 has ratio 1000/800=1.25 < 3, should be removed
        # Feature at 200.0 not in blank, should be kept
        mzs = [float(f["mz"]) for f in cleaned]
        assert 200.0 in mzs
        assert 100.0 not in mzs

    def test_keep_high_fold_change(self):
        from blank_subtraction_tool import subtract_blanks

        sample = [
            {"mz": "100.0", "rt": "60.0", "intensity": "10000"},
        ]
        blank = [
            {"mz": "100.0", "rt": "60.0", "intensity": "100"},
        ]

        cleaned = subtract_blanks(sample, blank, fold_change=3.0)
        # Ratio = 100, well above threshold
        assert len(cleaned) == 1

    def test_empty_blank(self):
        from blank_subtraction_tool import subtract_blanks

        sample = [
            {"mz": "100.0", "rt": "60.0", "intensity": "1000"},
        ]
        cleaned = subtract_blanks(sample, [], fold_change=3.0)
        assert len(cleaned) == 1

    def test_rt_tolerance(self):
        from blank_subtraction_tool import subtract_blanks

        sample = [
            {"mz": "100.0", "rt": "60.0", "intensity": "500"},
        ]
        blank = [
            {"mz": "100.0", "rt": "200.0", "intensity": "500"},  # far RT
        ]

        cleaned = subtract_blanks(sample, blank, fold_change=3.0, rt_tolerance=10.0)
        # RT difference 140s > 10s tolerance, so not matched
        assert len(cleaned) == 1

    def test_mz_tolerance(self):
        from blank_subtraction_tool import subtract_blanks

        sample = [
            {"mz": "100.0", "rt": "60.0", "intensity": "500"},
        ]
        blank = [
            {"mz": "100.001", "rt": "60.0", "intensity": "500"},
        ]

        # 0.001 Da at 100 m/z = 10 ppm, with tolerance 5 ppm should not match
        cleaned = subtract_blanks(sample, blank, fold_change=3.0, mz_tolerance_ppm=5.0)
        assert len(cleaned) == 1
