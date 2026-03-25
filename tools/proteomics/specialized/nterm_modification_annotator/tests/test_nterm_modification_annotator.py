"""Tests for nterm_modification_annotator."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestNtermModificationAnnotator:
    def _create_fasta(self, tmpdir, proteins):
        """Helper to create a FASTA file."""
        import pyopenms as oms

        fasta_path = os.path.join(tmpdir, "reference.fasta")
        entries = []
        for acc, seq in proteins.items():
            entry = oms.FASTAEntry()
            entry.identifier = acc
            entry.sequence = seq
            entries.append(entry)
        fasta_file = oms.FASTAFile()
        fasta_file.store(fasta_path, entries)
        return fasta_path

    def test_get_clean_sequence(self):
        from nterm_modification_annotator import get_clean_sequence
        assert get_clean_sequence("PEPTIDEK") == "PEPTIDEK"

    def test_detect_nterm_modification_acetyl(self):
        from nterm_modification_annotator import detect_nterm_modification
        assert detect_nterm_modification(".(Acetyl)PEPTIDEK") == "acetylation"

    def test_detect_nterm_modification_none(self):
        from nterm_modification_annotator import detect_nterm_modification
        assert detect_nterm_modification("PEPTIDEK") == "none"

    def test_detect_nterm_tmt(self):
        from nterm_modification_annotator import detect_nterm_modification
        assert detect_nterm_modification("TMT6plex-PEPTIDEK") == "TMT-label"

    def test_find_peptide_start(self):
        from nterm_modification_annotator import find_peptide_start
        assert find_peptide_start("DEF", "ABCDEFGHIJ") == 3
        assert find_peptide_start("ABC", "ABCDEFGHIJ") == 0
        assert find_peptide_start("ZZZ", "ABCDEFGHIJ") == -1

    def test_classify_protein_nterm(self):
        from nterm_modification_annotator import classify_nterm_type
        protein_seq = "MACDEFGHIJ"
        result = classify_nterm_type(0, protein_seq)
        assert result == "protein_nterm"

    def test_classify_met_removal(self):
        from nterm_modification_annotator import classify_nterm_type
        protein_seq = "MACDEFGHIJ"
        result = classify_nterm_type(1, protein_seq)
        assert result == "met_removal"

    def test_classify_neo_nterm(self):
        from nterm_modification_annotator import classify_nterm_type
        protein_seq = "MACDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZABCDEFGHIJKLMNOPQRSTUVWXYZ"
        # Position 75 - beyond signal peptide range
        result = classify_nterm_type(75, protein_seq)
        assert result == "neo_nterm"

    def test_classify_signal_peptide_known(self):
        from nterm_modification_annotator import classify_nterm_type
        protein_seq = "M" + "A" * 30 + "DEFGHIJ"
        sp_sites = {"P1": 25}
        result = classify_nterm_type(25, protein_seq, sp_sites, "P1")
        assert result == "signal_peptide"

    def test_classify_signal_peptide_candidate(self):
        from nterm_modification_annotator import classify_nterm_type
        # Position 20, preceded by 'A' (small residue)
        protein_seq = "M" + "K" * 18 + "A" + "DEFGHIJKLMNOP"
        result = classify_nterm_type(20, protein_seq)
        assert result == "signal_peptide_candidate"

    def test_classify_unmapped(self):
        from nterm_modification_annotator import classify_nterm_type
        result = classify_nterm_type(-1, "MACDEF")
        assert result == "unmapped"

    def test_annotate_nterm_peptides(self):
        from nterm_modification_annotator import annotate_nterm_peptides
        proteins = {"P1": "MACDEFGHIKLMNPQRSTVWY"}
        rows = [
            {"peptide": "MACDEF", "protein": "P1"},
            {"peptide": "ACDEF", "protein": "P1"},
            {"peptide": "GHIKLM", "protein": "P1"},
        ]
        results = annotate_nterm_peptides(rows, proteins)
        assert results[0]["nterm_type"] == "protein_nterm"
        assert results[1]["nterm_type"] == "met_removal"
        assert results[2]["nterm_type"] == "neo_nterm"

    def test_annotate_missing_protein(self):
        from nterm_modification_annotator import annotate_nterm_peptides
        rows = [{"peptide": "PEPTIDEK", "protein": "MISSING"}]
        results = annotate_nterm_peptides(rows, {})
        assert results[0]["nterm_type"] == "unmapped"

    def test_compute_summary(self):
        from nterm_modification_annotator import compute_summary
        results = [
            {"nterm_type": "protein_nterm"},
            {"nterm_type": "protein_nterm"},
            {"nterm_type": "neo_nterm"},
            {"nterm_type": "met_removal"},
        ]
        summary = compute_summary(results)
        assert summary["protein_nterm"] == 2
        assert summary["neo_nterm"] == 1
        assert summary["met_removal"] == 1

    def test_full_pipeline(self):
        from nterm_modification_annotator import annotate_nterm_peptides, load_fasta, write_output

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = self._create_fasta(tmpdir, {"P1": "MACDEFGHIKLMNPQRSTVWY"})
            proteins = load_fasta(fasta_path)
            rows = [
                {"peptide": "MACDEF", "protein": "P1"},
                {"peptide": "ACDEF", "protein": "P1"},
            ]
            results = annotate_nterm_peptides(rows, proteins)
            output_path = os.path.join(tmpdir, "annotated.tsv")
            write_output(output_path, results)
            assert os.path.exists(output_path)
            assert len(results) == 2
