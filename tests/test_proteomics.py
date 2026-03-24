"""
Tests for proteomics scripts.
Run with:  python -m pytest tests/ -v
"""

import sys
import os

import pytest

# Allow importing from the scripts directory
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts", "proteomics"))

try:
    import pyopenms as oms  # noqa: F401

    HAS_PYOPENMS = True
except ImportError:
    HAS_PYOPENMS = False

requires_pyopenms = pytest.mark.skipif(
    not HAS_PYOPENMS, reason="pyopenms not installed"
)


@requires_pyopenms
class TestPeptideMassCalculator:
    def test_basic_mass(self):
        from peptide_mass_calculator import peptide_masses

        result = peptide_masses("PEPTIDEK")
        assert result["sequence"] == "PEPTIDEK"
        assert result["charge"] == 1
        assert 927.0 < result["monoisotopic_mass"] < 928.0
        assert result["mz_monoisotopic"] > result["monoisotopic_mass"]

    def test_charge_state(self):
        from peptide_mass_calculator import peptide_masses

        r1 = peptide_masses("PEPTIDEK", charge=1)
        r2 = peptide_masses("PEPTIDEK", charge=2)
        # Higher charge → lower m/z
        assert r2["mz_monoisotopic"] < r1["mz_monoisotopic"]

    def test_fragment_ions(self):
        from peptide_mass_calculator import fragment_ions

        ions = fragment_ions("PEPTIDEK")
        seq_len = len("PEPTIDEK")
        # Expect seq_len-1 b-ions and seq_len-1 y-ions
        assert len(ions["b_ions"]) == seq_len - 1
        assert len(ions["y_ions"]) == seq_len - 1

    def test_modified_sequence(self):
        from peptide_mass_calculator import peptide_masses

        # Oxidised methionine
        result = peptide_masses("PEPTM[147]IDEK")
        assert result["monoisotopic_mass"] > 0

    def test_mz_formula(self):
        from peptide_mass_calculator import peptide_masses, PROTON

        r = peptide_masses("PEPTIDEK", charge=2)
        expected = (r["monoisotopic_mass"] + 2 * PROTON) / 2
        assert abs(r["mz_monoisotopic"] - expected) < 1e-6


@requires_pyopenms
class TestProteinDigest:
    # Short test protein containing known tryptic sites (K, R)
    PROTEIN = "MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELAL"

    def test_tryptic_digest_returns_peptides(self):
        from protein_digest import digest_protein

        peptides = digest_protein(self.PROTEIN, enzyme="Trypsin", min_length=1)
        assert len(peptides) > 0

    def test_peptide_structure(self):
        from protein_digest import digest_protein

        peptides = digest_protein(self.PROTEIN, enzyme="Trypsin", min_length=1)
        for pep in peptides:
            assert "sequence" in pep
            assert "monoisotopic_mass" in pep
            assert pep["monoisotopic_mass"] > 0

    def test_missed_cleavages(self):
        from protein_digest import digest_protein

        peps_0 = digest_protein(self.PROTEIN, enzyme="Trypsin", missed_cleavages=0, min_length=1)
        peps_2 = digest_protein(self.PROTEIN, enzyme="Trypsin", missed_cleavages=2, min_length=1)
        # More missed cleavages → at least as many peptides
        assert len(peps_2) >= len(peps_0)

    def test_length_filter(self):
        from protein_digest import digest_protein

        peptides = digest_protein(
            self.PROTEIN, enzyme="Trypsin", min_length=5, max_length=20, missed_cleavages=2
        )
        for pep in peptides:
            assert 5 <= pep["length"] <= 20

    def test_list_enzymes(self):
        from protein_digest import list_enzymes

        enzymes = list_enzymes()
        assert "Trypsin" in enzymes
        assert len(enzymes) > 5
