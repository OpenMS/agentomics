"""Tests for duplicate_feature_detector."""

import pytest

pytest.importorskip("pyopenms")


class TestDuplicateFeatureDetector:
    def test_detect_duplicates(self):
        from duplicate_feature_detector import detect_duplicates

        features = [
            {"mz": "100.0000", "rt": "60.0", "intensity": "1000"},
            {"mz": "100.0005", "rt": "61.0", "intensity": "800"},  # duplicate
            {"mz": "200.0000", "rt": "120.0", "intensity": "500"},
        ]
        result = detect_duplicates(features, mz_tolerance_ppm=10.0, rt_tolerance=5.0)
        assert len(result) == 3
        # First two should share a group
        assert result[0]["group_id"] == result[1]["group_id"]
        assert result[0]["group_id"] != result[2]["group_id"]

    def test_deduplicate(self):
        from duplicate_feature_detector import deduplicate

        features = [
            {"mz": "100.0000", "rt": "60.0", "intensity": "1000"},
            {"mz": "100.0005", "rt": "61.0", "intensity": "800"},
            {"mz": "200.0000", "rt": "120.0", "intensity": "500"},
        ]
        deduped = deduplicate(features, mz_tolerance_ppm=10.0, rt_tolerance=5.0)
        assert len(deduped) == 2

    def test_keep_highest_intensity(self):
        from duplicate_feature_detector import deduplicate

        features = [
            {"mz": "100.0000", "rt": "60.0", "intensity": "500"},
            {"mz": "100.0005", "rt": "61.0", "intensity": "1000"},  # higher intensity
        ]
        deduped = deduplicate(features, mz_tolerance_ppm=10.0, rt_tolerance=5.0)
        assert len(deduped) == 1
        assert float(deduped[0]["intensity"]) == 1000.0

    def test_no_duplicates(self):
        from duplicate_feature_detector import deduplicate

        features = [
            {"mz": "100.0", "rt": "60.0", "intensity": "1000"},
            {"mz": "500.0", "rt": "300.0", "intensity": "500"},
        ]
        deduped = deduplicate(features, mz_tolerance_ppm=10.0, rt_tolerance=5.0)
        assert len(deduped) == 2

    def test_all_duplicates(self):
        from duplicate_feature_detector import deduplicate

        features = [
            {"mz": "100.0000", "rt": "60.0", "intensity": "1000"},
            {"mz": "100.0001", "rt": "60.1", "intensity": "999"},
            {"mz": "100.0002", "rt": "60.2", "intensity": "998"},
        ]
        deduped = deduplicate(features, mz_tolerance_ppm=10.0, rt_tolerance=5.0)
        assert len(deduped) == 1
