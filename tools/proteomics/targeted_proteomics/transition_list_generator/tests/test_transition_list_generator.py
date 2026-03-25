"""Tests for transition_list_generator."""

import pytest

pytest.importorskip("pyopenms")


class TestTransitionListGenerator:
    def test_basic_transitions(self):
        from transition_list_generator import generate_transitions

        transitions = generate_transitions("PEPTIDEK", precursor_charges=[2])
        assert len(transitions) > 0
        for t in transitions:
            assert t["precursor_mz"] > 0
            assert t["product_mz"] > 0
            assert t["precursor_charge"] == 2

    def test_multiple_charges(self):
        from transition_list_generator import generate_transitions

        transitions = generate_transitions("PEPTIDEK", precursor_charges=[2, 3])
        charges = set(t["precursor_charge"] for t in transitions)
        assert 2 in charges
        assert 3 in charges

    def test_ion_range_filter(self):
        from transition_list_generator import generate_transitions

        all_trans = generate_transitions("PEPTIDEK", precursor_charges=[2])
        filtered = generate_transitions("PEPTIDEK", precursor_charges=[2], ion_range="y3-y6")
        assert len(filtered) <= len(all_trans)
        for t in filtered:
            assert t["annotation"].startswith("y")

    def test_parse_ion_range(self):
        from transition_list_generator import parse_ion_range

        series, start, end = parse_ion_range("y3-y8")
        assert series == "y"
        assert start == 3
        assert end == 8

    def test_precursor_mz_decreases_with_charge(self):
        from transition_list_generator import generate_transitions

        t2 = generate_transitions("PEPTIDEK", precursor_charges=[2])
        t3 = generate_transitions("PEPTIDEK", precursor_charges=[3])
        assert t2[0]["precursor_mz"] > t3[0]["precursor_mz"]

    def test_transitions_have_annotations(self):
        from transition_list_generator import generate_transitions

        transitions = generate_transitions("PEPTIDEK", precursor_charges=[2])
        annotated = [t for t in transitions if t["annotation"]]
        assert len(annotated) > 0
