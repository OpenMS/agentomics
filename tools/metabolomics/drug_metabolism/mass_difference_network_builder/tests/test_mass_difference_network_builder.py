"""Tests for mass_difference_network_builder."""

import pytest

pytest.importorskip("pyopenms")


class TestMassDifferenceNetworkBuilder:
    def test_oxidation_edge(self):
        from mass_difference_network_builder import build_network

        features = [
            {"feature_id": "A", "mz": "200.0000"},
            {"feature_id": "B", "mz": "215.9949"},  # +15.9949 = oxidation
        ]
        reactions = [("Oxidation", 15.994915)]
        edges = build_network(features, reactions, tolerance=0.005)
        assert len(edges) >= 1
        rxn_names = [e["reaction"] for e in edges]
        assert "Oxidation" in rxn_names

    def test_no_match(self):
        from mass_difference_network_builder import build_network

        features = [
            {"feature_id": "A", "mz": "200.0"},
            {"feature_id": "B", "mz": "300.0"},
        ]
        reactions = [("Oxidation", 15.994915)]
        edges = build_network(features, reactions, tolerance=0.005)
        assert len(edges) == 0

    def test_multiple_reactions(self):
        from mass_difference_network_builder import build_network

        features = [
            {"feature_id": "A", "mz": "200.0000"},
            {"feature_id": "B", "mz": "215.9949"},
            {"feature_id": "C", "mz": "214.0157"},  # +14.0157 = methylation
        ]
        reactions = [("Oxidation", 15.994915), ("Methylation", 14.015650)]
        edges = build_network(features, reactions, tolerance=0.005)
        rxn_names = set(e["reaction"] for e in edges)
        assert "Oxidation" in rxn_names
        assert "Methylation" in rxn_names

    def test_default_reactions(self):
        from mass_difference_network_builder import DEFAULT_REACTIONS

        assert len(DEFAULT_REACTIONS) > 10
        names = [r[0] for r in DEFAULT_REACTIONS]
        assert "Oxidation" in names
        assert "Dehydration" in names

    def test_load_reactions_none(self):
        from mass_difference_network_builder import DEFAULT_REACTIONS, load_reactions

        reactions = load_reactions(None)
        assert reactions == DEFAULT_REACTIONS
