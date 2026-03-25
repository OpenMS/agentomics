"""Tests for adduct_group_analyzer."""

import pytest

pytest.importorskip("pyopenms")


class TestAdductGroupAnalyzer:
    def test_group_mh_mna(self):
        from adduct_group_analyzer import PROTON, find_adduct_groups

        # [M+H]+ and [M+Na]+ differ by 22.989218 - 1.007276 = 21.981942
        mh_mz = 181.0707  # glucose [M+H]+
        mna_mz = mh_mz + (22.989218 - PROTON)

        features = [
            {"feature_id": "0", "mz": str(mh_mz), "rt": "100.0"},
            {"feature_id": "1", "mz": str(mna_mz), "rt": "100.0"},
            {"feature_id": "2", "mz": "500.0", "rt": "200.0"},
        ]

        groups = find_adduct_groups(features, rt_tolerance=5.0, mz_tolerance_da=0.02)
        # Features 0 and 1 should share a group
        g0 = next(g for g in groups if g["feature_id"] == "0")["group_id"]
        g1 = next(g for g in groups if g["feature_id"] == "1")["group_id"]
        g2 = next(g for g in groups if g["feature_id"] == "2")["group_id"]
        assert g0 == g1
        assert g0 != g2

    def test_rt_separation_prevents_grouping(self):
        from adduct_group_analyzer import PROTON, find_adduct_groups

        mh_mz = 181.0707
        mna_mz = mh_mz + (22.989218 - PROTON)

        features = [
            {"feature_id": "0", "mz": str(mh_mz), "rt": "100.0"},
            {"feature_id": "1", "mz": str(mna_mz), "rt": "200.0"},  # far RT
        ]

        groups = find_adduct_groups(features, rt_tolerance=5.0)
        g0 = next(g for g in groups if g["feature_id"] == "0")["group_id"]
        g1 = next(g for g in groups if g["feature_id"] == "1")["group_id"]
        assert g0 != g1

    def test_single_feature(self):
        from adduct_group_analyzer import find_adduct_groups

        features = [{"feature_id": "0", "mz": "500.0", "rt": "100.0"}]
        groups = find_adduct_groups(features)
        assert len(groups) == 1

    def test_all_features_annotated(self):
        from adduct_group_analyzer import find_adduct_groups

        features = [
            {"feature_id": str(i), "mz": str(200.0 + i * 50), "rt": "100.0"}
            for i in range(5)
        ]
        groups = find_adduct_groups(features)
        assert len(groups) == 5
        assert all("group_id" in g for g in groups)
