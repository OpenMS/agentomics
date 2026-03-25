"""Tests for peptide_uniqueness_checker."""

import tempfile

from conftest import requires_pyopenms


@requires_pyopenms
class TestPeptideUniquenessChecker:
    def _create_fasta(self, tmpdir):
        """Create a synthetic FASTA file with two proteins."""
        import pyopenms as oms

        fasta_path = f"{tmpdir}/test.fasta"
        entries = []
        e1 = oms.FASTAEntry()
        e1.identifier = "sp|P00001|PROT1"
        e1.sequence = "MSPEPTIDEKAAANOTHERPEPTIDE"
        entries.append(e1)
        e2 = oms.FASTAEntry()
        e2.identifier = "sp|P00002|PROT2"
        e2.sequence = "MSPEPTIDEKGGGUNIQUEPEPTIDE"
        entries.append(e2)
        oms.FASTAFile().store(fasta_path, entries)
        return fasta_path

    def test_proteotypic_peptide(self):
        from peptide_uniqueness_checker import check_uniqueness

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = self._create_fasta(tmpdir)
            results = check_uniqueness(["UNIQUEPEPTIDE"], fasta_path)
            assert len(results) == 1
            assert results[0]["is_proteotypic"] is True
            assert results[0]["protein_count"] == 1

    def test_shared_peptide(self):
        from peptide_uniqueness_checker import check_uniqueness

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = self._create_fasta(tmpdir)
            results = check_uniqueness(["PEPTIDEK"], fasta_path)
            assert len(results) == 1
            assert results[0]["is_proteotypic"] is False
            assert results[0]["protein_count"] == 2

    def test_missing_peptide(self):
        from peptide_uniqueness_checker import check_uniqueness

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = self._create_fasta(tmpdir)
            results = check_uniqueness(["ZZZZZZZ"], fasta_path)
            assert results[0]["protein_count"] == 0
            assert results[0]["is_proteotypic"] is False

    def test_load_fasta(self):
        from peptide_uniqueness_checker import load_fasta

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = self._create_fasta(tmpdir)
            entries = load_fasta(fasta_path)
            assert len(entries) == 2
