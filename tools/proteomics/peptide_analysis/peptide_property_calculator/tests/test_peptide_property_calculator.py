"""Tests for peptide_property_calculator."""

import pytest

pytest.importorskip("pyopenms")


class TestPeptidePropertyCalculator:
    def test_calculate_properties_basic(self):
        from peptide_property_calculator import calculate_properties

        result = calculate_properties("PEPTIDEK")
        assert result["unmodified_sequence"] == "PEPTIDEK"
        assert result["length"] == 8
        assert result["monoisotopic_mass"] > 900
        assert 3.0 < result["pI"] < 7.0
        assert isinstance(result["gravy"], float)
        assert isinstance(result["charge_at_ph"], float)

    def test_pi_basic_peptide(self):
        from peptide_property_calculator import calculate_pi

        pi = calculate_pi("ACDEFGHIK")
        assert 4.0 < pi < 7.0

    def test_gravy(self):
        from peptide_property_calculator import calculate_gravy

        gravy = calculate_gravy("AAAA")
        assert gravy == 1.8  # all alanine

    def test_charge_at_ph(self):
        from peptide_property_calculator import _charge_at_ph

        # At very low pH, charge should be positive
        assert _charge_at_ph("PEPTIDEK", 1.0) > 0
        # At very high pH, charge should be negative
        assert _charge_at_ph("PEPTIDEK", 14.0) < 0

    def test_amino_acid_composition(self):
        from peptide_property_calculator import amino_acid_composition

        comp = amino_acid_composition("AAAK")
        assert comp["counts"]["A"] == 3
        assert comp["counts"]["K"] == 1
        assert abs(comp["frequencies"]["A"] - 0.75) < 0.01

    def test_instability_index(self):
        from peptide_property_calculator import calculate_instability_index

        ii = calculate_instability_index("DGDG")
        assert ii > 0

    def test_output_json(self):
        import json
        import tempfile

        from peptide_property_calculator import calculate_properties

        result = calculate_properties("PEPTIDEK")
        with tempfile.NamedTemporaryFile(suffix=".json", mode="w", delete=False) as fh:
            json.dump(result, fh)
