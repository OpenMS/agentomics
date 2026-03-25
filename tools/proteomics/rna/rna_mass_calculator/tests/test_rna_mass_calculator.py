"""Tests for rna_mass_calculator."""


import pytest
from conftest import requires_pyopenms
from rna_mass_calculator import _manual_formula, calculate_isotope_pattern, calculate_rna_mass


@requires_pyopenms
class TestRnaMassCalculator:
    def test_basic_mass(self):
        result = calculate_rna_mass("AAUGC", charge=1)
        assert result["sequence"] == "AAUGC"
        assert result["charge"] == 1
        assert result["monoisotopic_mass"] > 0
        assert result["mz"] > 0
        assert len(result["formula"]) > 0

    def test_charge_state(self):
        r1 = calculate_rna_mass("AAUGC", charge=1)
        r2 = calculate_rna_mass("AAUGC", charge=2)
        assert r1["monoisotopic_mass"] == r2["monoisotopic_mass"]
        assert r2["mz"] < r1["mz"]

    def test_longer_sequence(self):
        result = calculate_rna_mass("AAUGCAAUGG", charge=3)
        assert result["charge"] == 3
        assert result["monoisotopic_mass"] > 0

    def test_invalid_nucleotide(self):
        with pytest.raises(ValueError, match="Invalid RNA nucleotide"):
            calculate_rna_mass("AATGC")

    def test_case_insensitive(self):
        r1 = calculate_rna_mass("aaugc")
        r2 = calculate_rna_mass("AAUGC")
        assert abs(r1["monoisotopic_mass"] - r2["monoisotopic_mass"]) < 0.001

    def test_isotope_pattern(self):
        pattern = calculate_isotope_pattern("AAUGC", n_peaks=5)
        assert len(pattern) == 5
        assert all(m > 0 for m, _ in pattern)
        assert all(0 <= i <= 1.0 for _, i in pattern)

    def test_manual_formula(self):
        formula = _manual_formula("AU")
        assert "C" in formula
        assert "H" in formula
        assert "N" in formula
        assert "O" in formula
        assert "P" in formula
