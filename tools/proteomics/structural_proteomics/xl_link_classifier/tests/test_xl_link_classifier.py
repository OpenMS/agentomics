"""Tests for xl_link_classifier."""

import os
import tempfile

from conftest import requires_pyopenms


@requires_pyopenms
class TestXlLinkClassifier:
    def _create_fasta(self, tmpdir, proteins):
        """Helper to create a FASTA file."""
        import pyopenms as oms

        fasta_path = os.path.join(tmpdir, "proteome.fasta")
        entries = []
        for acc, seq in proteins.items():
            entry = oms.FASTAEntry()
            entry.identifier = acc
            entry.sequence = seq
            entries.append(entry)
        fasta_file = oms.FASTAFile()
        fasta_file.store(fasta_path, entries)
        return fasta_path

    def test_load_fasta(self):
        from xl_link_classifier import load_fasta
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = self._create_fasta(tmpdir, {"P1": "ACDEFGHIKLMNPQRSTVWY"})
            proteins = load_fasta(fasta_path)
            assert "P1" in proteins

    def test_strip_modifications(self):
        from xl_link_classifier import strip_modifications
        assert strip_modifications("PEPTM[147]IDEK") == "PEPTMIDEK"
        assert strip_modifications("PEPT(Oxidation)IDEK") == "PEPTIDEK"

    def test_find_peptide_proteins(self):
        from xl_link_classifier import find_peptide_proteins
        proteins = {"P1": "ACDEFGHIKLMNPQRSTVWY", "P2": "XXACDEFYYYYY"}
        found = find_peptide_proteins("ACDEF", proteins)
        assert "P1" in found
        assert "P2" in found

    def test_find_peptide_not_found(self):
        from xl_link_classifier import find_peptide_proteins
        proteins = {"P1": "ACDEFGHIK"}
        found = find_peptide_proteins("ZZZZZ", proteins)
        assert len(found) == 0

    def test_classify_monolink(self):
        from xl_link_classifier import classify_crosslink
        proteins = {"P1": "ACDEFGHIKLMNPQRSTVWY"}
        result = classify_crosslink("ACDEF", "", proteins)
        assert result["link_type"] == "monolink"

    def test_classify_monolink_dash(self):
        from xl_link_classifier import classify_crosslink
        proteins = {"P1": "ACDEFGHIKLMNPQRSTVWY"}
        result = classify_crosslink("ACDEF", "-", proteins)
        assert result["link_type"] == "monolink"

    def test_classify_intra_protein(self):
        from xl_link_classifier import classify_crosslink
        proteins = {"P1": "ACDEFGHIKLMNPQRSTVWY"}
        result = classify_crosslink("ACDEF", "GHIKLM", proteins)
        assert result["link_type"] == "intra-protein"

    def test_classify_inter_protein(self):
        from xl_link_classifier import classify_crosslink
        proteins = {"P1": "ACDEFGHIK", "P2": "LMNPQRSTV"}
        result = classify_crosslink("ACDEF", "LMNPQ", proteins)
        assert result["link_type"] == "inter-protein"

    def test_classify_unknown_unmapped(self):
        from xl_link_classifier import classify_crosslink
        proteins = {"P1": "ACDEFGHIK"}
        result = classify_crosslink("ZZZZZ", "YYYYY", proteins)
        assert result["link_type"] == "unknown"

    def test_classify_crosslinks_batch(self):
        from xl_link_classifier import classify_crosslinks
        proteins = {"P1": "ACDEFGHIKLMNPQRSTVWY", "P2": "XXXXXXYYYYYY"}
        crosslinks = [
            {"peptide1": "ACDEF", "peptide2": "GHIKLM"},
            {"peptide1": "ACDEF", "peptide2": "XXXXXX"},
            {"peptide1": "ACDEF", "peptide2": ""},
        ]
        results = classify_crosslinks(crosslinks, proteins)
        assert results[0]["link_type"] == "intra-protein"
        assert results[1]["link_type"] == "inter-protein"
        assert results[2]["link_type"] == "monolink"

    def test_compute_summary(self):
        from xl_link_classifier import compute_summary
        results = [
            {"link_type": "intra-protein"},
            {"link_type": "intra-protein"},
            {"link_type": "inter-protein"},
            {"link_type": "monolink"},
        ]
        summary = compute_summary(results)
        assert summary["intra-protein"] == 2
        assert summary["inter-protein"] == 1
        assert summary["monolink"] == 1

    def test_full_pipeline(self):
        from xl_link_classifier import classify_crosslinks, load_fasta, write_output

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = self._create_fasta(tmpdir, {
                "P1": "ACDEFGHIKLMNPQRSTVWY",
                "P2": "XXXXXXYYYYYY",
            })
            proteins = load_fasta(fasta_path)
            crosslinks = [
                {"peptide1": "ACDEF", "peptide2": "GHIKLM"},
                {"peptide1": "ACDEF", "peptide2": "XXXXXX"},
            ]
            results = classify_crosslinks(crosslinks, proteins)
            output_path = os.path.join(tmpdir, "classified.tsv")
            write_output(output_path, results)
            assert os.path.exists(output_path)
