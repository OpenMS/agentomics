"""Tests for protein_coverage_calculator."""

import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestProteinCoverageCalculator:
    def _create_fasta(self, tmpdir):
        import pyopenms as oms

        fasta_path = f"{tmpdir}/test.fasta"
        entries = []
        e1 = oms.FASTAEntry()
        e1.identifier = "PROT1"
        e1.sequence = "MSPEPTIDEKAAANOTHERPEPTIDE"
        entries.append(e1)
        oms.FASTAFile().store(fasta_path, entries)
        return fasta_path

    def test_full_coverage(self):
        from protein_coverage_calculator import calculate_coverage

        proteins = {"PROT1": "PEPTIDEK"}
        results = calculate_coverage(proteins, ["PEPTIDEK"])
        assert results[0]["coverage_percent"] == 100.0

    def test_partial_coverage(self):
        from protein_coverage_calculator import calculate_coverage

        proteins = {"PROT1": "MSPEPTIDEKAAAA"}
        results = calculate_coverage(proteins, ["PEPTIDEK"])
        assert 0 < results[0]["coverage_percent"] < 100

    def test_no_match(self):
        from protein_coverage_calculator import calculate_coverage

        proteins = {"PROT1": "MSPEPTIDEK"}
        results = calculate_coverage(proteins, ["ZZZZZ"])
        assert results[0]["coverage_percent"] == 0.0

    def test_load_fasta(self):
        from protein_coverage_calculator import load_fasta

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = self._create_fasta(tmpdir)
            proteins = load_fasta(fasta_path)
            assert "PROT1" in proteins

    def test_multiple_peptides(self):
        from protein_coverage_calculator import calculate_coverage

        proteins = {"PROT1": "MSPEPTIDEKAAANOTHERPEPTIDE"}
        results = calculate_coverage(proteins, ["PEPTIDEK", "ANOTHERPEPTIDE"])
        assert results[0]["matched_peptides"] == 2
        assert results[0]["coverage_percent"] > 50
