"""Tests for isf_detector."""

import pytest
from conftest import requires_pyopenms


@requires_pyopenms
class TestNeutralLossMasses:
    def test_h2o_mass(self):
        from isf_detector import get_neutral_loss_masses

        masses = get_neutral_loss_masses()
        assert masses["H2O"] == pytest.approx(18.0106, abs=0.001)

    def test_co2_mass(self):
        from isf_detector import get_neutral_loss_masses

        masses = get_neutral_loss_masses()
        assert masses["CO2"] == pytest.approx(43.9898, abs=0.001)


@requires_pyopenms
class TestPearsonCorrelation:
    def test_perfect_correlation(self):
        from isf_detector import pearson_correlation

        assert pearson_correlation([1, 2, 3], [2, 4, 6]) == pytest.approx(1.0, abs=0.001)

    def test_no_correlation(self):
        from isf_detector import pearson_correlation

        r = pearson_correlation([1, 2, 3, 4], [1, -1, 1, -1])
        assert abs(r) < 0.5

    def test_short_input(self):
        from isf_detector import pearson_correlation

        assert pearson_correlation([1], [2]) == 0.0


@requires_pyopenms
class TestDetectISFPairs:
    def test_water_loss_detected(self):
        from isf_detector import detect_isf_pairs, get_neutral_loss_masses

        h2o_mass = get_neutral_loss_masses()["H2O"]
        features = [
            {"id": "F1", "mz": "180.0634", "rt": "120.5", "intensity": "50000"},
            {"id": "F2", "mz": str(180.0634 - h2o_mass), "rt": "121.0", "intensity": "15000"},
        ]
        pairs = detect_isf_pairs(features, rt_tolerance=3.0, mass_tolerance_da=0.01)
        assert len(pairs) == 1
        assert pairs[0]["neutral_loss"] == "H2O"
        assert pairs[0]["precursor_id"] == "F1"

    def test_no_pair_when_rt_too_far(self):
        from isf_detector import detect_isf_pairs, get_neutral_loss_masses

        h2o_mass = get_neutral_loss_masses()["H2O"]
        features = [
            {"id": "F1", "mz": "180.0634", "rt": "100.0", "intensity": "50000"},
            {"id": "F2", "mz": str(180.0634 - h2o_mass), "rt": "200.0", "intensity": "15000"},
        ]
        pairs = detect_isf_pairs(features, rt_tolerance=3.0)
        assert len(pairs) == 0

    def test_no_pair_when_mass_not_matching(self):
        from isf_detector import detect_isf_pairs

        features = [
            {"id": "F1", "mz": "180.0634", "rt": "120.5", "intensity": "50000"},
            {"id": "F2", "mz": "170.0000", "rt": "121.0", "intensity": "15000"},
        ]
        pairs = detect_isf_pairs(features, rt_tolerance=3.0)
        assert len(pairs) == 0


@requires_pyopenms
class TestAnnotateFeatures:
    def test_annotation(self):
        from isf_detector import annotate_features

        features = [
            {"id": "F1", "mz": "180.0", "rt": "120.5", "intensity": "50000"},
            {"id": "F2", "mz": "162.0", "rt": "121.0", "intensity": "15000"},
            {"id": "F3", "mz": "300.0", "rt": "200.0", "intensity": "80000"},
        ]
        pairs = [{"precursor_id": "F1", "fragment_id": "F2", "neutral_loss": "H2O"}]
        annotated = annotate_features(features, pairs)

        assert annotated[0]["isf_flag"] is True
        assert annotated[0]["isf_role"] == "precursor"
        assert annotated[1]["isf_flag"] is True
        assert annotated[1]["isf_role"] == "fragment"
        assert annotated[2]["isf_flag"] is False
