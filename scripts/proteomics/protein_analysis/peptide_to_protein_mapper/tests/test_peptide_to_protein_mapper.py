"""Tests for peptide_to_protein_mapper."""


from conftest import requires_pyopenms
from peptide_to_protein_mapper import _strip_modifications, map_peptides_to_proteins, read_fasta


@requires_pyopenms
class TestPeptideToProteinMapper:
    def _make_fasta_entries(self):
        return [
            ("sp|P12345|PROT1", "Protein 1", "MKPEPTIDEKHELLO"),
            ("sp|P67890|PROT2", "Protein 2", "TESTPEPWORLDANOTHERPEP"),
            ("sp|P11111|PROT3", "Protein 3", "MKPEPTIDEKWORLD"),  # shares PEPTIDEK
        ]

    def test_basic_mapping(self):
        fasta = self._make_fasta_entries()
        peptides = [{"peptide": "PEPTIDEK"}]
        result = map_peptides_to_proteins(peptides, fasta)
        proteins = [r["protein"] for r in result]
        assert "sp|P12345|PROT1" in proteins
        assert "sp|P11111|PROT3" in proteins

    def test_unique_peptide(self):
        fasta = self._make_fasta_entries()
        peptides = [{"peptide": "TESTPEP"}]
        result = map_peptides_to_proteins(peptides, fasta)
        assert len(result) == 1
        assert result[0]["is_unique"] is True

    def test_shared_peptide(self):
        fasta = self._make_fasta_entries()
        peptides = [{"peptide": "PEPTIDEK"}]
        result = map_peptides_to_proteins(peptides, fasta)
        assert len(result) == 2
        assert all(not r["is_unique"] for r in result)

    def test_unmapped_peptide(self):
        fasta = self._make_fasta_entries()
        peptides = [{"peptide": "NONEXISTENT"}]
        result = map_peptides_to_proteins(peptides, fasta)
        assert len(result) == 1
        assert result[0]["protein"] == ""

    def test_start_end_positions(self):
        fasta = self._make_fasta_entries()
        peptides = [{"peptide": "PEPTIDEK"}]
        result = map_peptides_to_proteins(peptides, fasta)
        p1_mapping = next(r for r in result if r["protein"] == "sp|P12345|PROT1")
        assert p1_mapping["start"] == 3  # 1-based, MK|PEPTIDEK|HELLO
        assert p1_mapping["end"] == 10

    def test_strip_modifications(self):
        assert _strip_modifications("PEPTM[147]IDEK") == "PEPTMIDEK"
        assert _strip_modifications("PEPTIDEK") == "PEPTIDEK"
        assert _strip_modifications("[Acetyl]PEPTIDEK") == "PEPTIDEK"

    def test_read_fasta(self, tmp_path):
        fasta_file = str(tmp_path / "test.fasta")
        with open(fasta_file, "w") as fh:
            fh.write(">sp|P12345|PROT1 Test protein\n")
            fh.write("MKPEPTIDEK\n")
        entries = read_fasta(fasta_file)
        assert len(entries) == 1
        assert "P12345" in entries[0][0]
