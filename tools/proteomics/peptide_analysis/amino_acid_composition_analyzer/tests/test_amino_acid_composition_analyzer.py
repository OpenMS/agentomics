"""Tests for amino_acid_composition_analyzer."""

import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestAminoAcidCompositionAnalyzer:
    def test_basic_composition(self):
        from amino_acid_composition_analyzer import analyze_composition

        result = analyze_composition("AAAKKK")
        assert result["length"] == 6
        assert result["counts"]["A"] == 3
        assert result["counts"]["K"] == 3
        assert abs(result["frequencies"]["A"] - 0.5) < 0.01

    def test_all_standard_aas(self):
        from amino_acid_composition_analyzer import STANDARD_AAS, analyze_composition

        result = analyze_composition("ACDEFGHIKLMNPQRSTVWY")
        for aa in STANDARD_AAS:
            assert result["counts"][aa] == 1

    def test_group_counts(self):
        from amino_acid_composition_analyzer import analyze_composition

        result = analyze_composition("KRHDE")
        assert result["basic_residues"] == 3  # K, R, H
        assert result["acidic_residues"] == 2  # D, E

    def test_hydrophobic_count(self):
        from amino_acid_composition_analyzer import analyze_composition

        result = analyze_composition("AILMFWV")
        assert result["hydrophobic_residues"] == 7

    def test_mass_positive(self):
        from amino_acid_composition_analyzer import analyze_composition

        result = analyze_composition("PEPTIDEK")
        assert result["monoisotopic_mass"] > 0

    def test_fasta_analysis(self):
        import pyopenms as oms
        from amino_acid_composition_analyzer import analyze_fasta

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = f"{tmpdir}/test.fasta"
            entries = []
            e1 = oms.FASTAEntry()
            e1.identifier = "PROT1"
            e1.sequence = "MSPEPTIDEK"
            entries.append(e1)
            e2 = oms.FASTAEntry()
            e2.identifier = "PROT2"
            e2.sequence = "ANOTHERPEPTIDE"
            entries.append(e2)
            oms.FASTAFile().store(fasta_path, entries)

            results = analyze_fasta(fasta_path)
            assert len(results) == 2
            assert results[0]["accession"] == "PROT1"
            assert results[1]["accession"] == "PROT2"

    def test_empty_counts_for_missing_aa(self):
        from amino_acid_composition_analyzer import analyze_composition

        result = analyze_composition("AAAA")
        assert result["counts"]["W"] == 0
        assert result["frequencies"]["W"] == 0.0
