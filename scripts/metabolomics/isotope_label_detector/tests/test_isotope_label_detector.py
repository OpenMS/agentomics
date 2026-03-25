"""Tests for isotope_label_detector."""

import pytest
from conftest import requires_pyopenms


@requires_pyopenms
class TestGetElementCount:
    def test_glucose_carbons(self):
        from isotope_label_detector import get_element_count

        assert get_element_count("C6H12O6", "C") == 6

    def test_alanine_nitrogen(self):
        from isotope_label_detector import get_element_count

        assert get_element_count("C3H7NO2", "N") == 1


@requires_pyopenms
class TestComputeExpectedMassShift:
    def test_glucose_13c(self):
        from isotope_label_detector import compute_expected_mass_shift

        shift = compute_expected_mass_shift("C6H12O6", "13C")
        # 6 carbons * 1.003355 Da
        assert shift == pytest.approx(6 * 1.003355, abs=0.001)

    def test_unsupported_tracer(self):
        from isotope_label_detector import compute_expected_mass_shift

        with pytest.raises(ValueError, match="Unsupported tracer"):
            compute_expected_mass_shift("C6H12O6", "18O")


@requires_pyopenms
class TestFindLabeledPairs:
    def test_match_13c_glucose(self):
        from isotope_label_detector import find_labeled_pairs

        shift_6c = 6 * 1.003355
        unlabeled = [{"id": "U1", "mz": "180.0634", "rt": "120.0"}]
        labeled = [{"id": "L1", "mz": str(180.0634 + shift_6c), "rt": "120.5"}]

        pairs = find_labeled_pairs(unlabeled, labeled, tracer="13C", ppm=5.0, rt_tolerance=10.0)
        assert len(pairs) == 1
        assert pairs[0]["n_labels"] == 6
        assert pairs[0]["ppm_error"] < 5.0

    def test_no_match_wrong_mass(self):
        from isotope_label_detector import find_labeled_pairs

        unlabeled = [{"id": "U1", "mz": "180.0634", "rt": "120.0"}]
        labeled = [{"id": "L1", "mz": "200.0000", "rt": "120.5"}]

        pairs = find_labeled_pairs(unlabeled, labeled, tracer="13C", ppm=5.0)
        assert len(pairs) == 0

    def test_no_match_rt_too_far(self):
        from isotope_label_detector import find_labeled_pairs

        shift_6c = 6 * 1.003355
        unlabeled = [{"id": "U1", "mz": "180.0634", "rt": "100.0"}]
        labeled = [{"id": "L1", "mz": str(180.0634 + shift_6c), "rt": "200.0"}]

        pairs = find_labeled_pairs(unlabeled, labeled, tracer="13C", ppm=5.0, rt_tolerance=10.0)
        assert len(pairs) == 0

    def test_multiple_pairs(self):
        from isotope_label_detector import find_labeled_pairs

        shift_3c = 3 * 1.003355
        shift_6c = 6 * 1.003355
        unlabeled = [
            {"id": "U1", "mz": "89.0477", "rt": "50.0"},
            {"id": "U2", "mz": "180.0634", "rt": "120.0"},
        ]
        labeled = [
            {"id": "L1", "mz": str(89.0477 + shift_3c), "rt": "50.5"},
            {"id": "L2", "mz": str(180.0634 + shift_6c), "rt": "120.5"},
        ]

        pairs = find_labeled_pairs(unlabeled, labeled, tracer="13C", ppm=5.0, rt_tolerance=10.0)
        assert len(pairs) == 2
