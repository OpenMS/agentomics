"""Tests for protein_digest."""


import pytest

pytest.importorskip("pyopenms")


class TestProteinDigest:
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
