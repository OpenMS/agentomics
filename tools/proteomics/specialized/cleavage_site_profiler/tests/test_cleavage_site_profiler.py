"""Tests for cleavage_site_profiler."""

import os
import tempfile

from conftest import requires_pyopenms


@requires_pyopenms
class TestCleavageSiteProfiler:
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
        from cleavage_site_profiler import get_clean_sequence
        assert get_clean_sequence("PEPTIDEK") == "PEPTIDEK"

    def test_find_cleavage_position(self):
        from cleavage_site_profiler import find_cleavage_position
        protein_seq = "ABCDEFGHIJKLMNOP"
        # Peptide "GHIJ" starts at position 6 in the protein
        pos = find_cleavage_position("GHIJ", "P1", protein_seq)
        assert pos == 6

    def test_find_cleavage_position_at_start(self):
        from cleavage_site_profiler import find_cleavage_position
        protein_seq = "ABCDEFGHIJ"
        # Peptide at position 0 -> no cleavage site (it's the native N-term)
        pos = find_cleavage_position("ABCD", "P1", protein_seq)
        assert pos == -1

    def test_find_cleavage_position_not_found(self):
        from cleavage_site_profiler import find_cleavage_position
        pos = find_cleavage_position("ZZZZZ", "P1", "ABCDEFGH")
        assert pos == -1

    def test_extract_cleavage_window_center(self):
        from cleavage_site_profiler import extract_cleavage_window
        protein_seq = "ABCDEFGHIJKLMNOP"
        # Cleavage at position 8 (between H and I)
        window = extract_cleavage_window(protein_seq, 8, 4)
        assert window == "EFGHIJKL"
        assert len(window) == 8

    def test_extract_cleavage_window_near_start(self):
        from cleavage_site_profiler import extract_cleavage_window
        protein_seq = "ABCDEFGHIJ"
        # Cleavage at position 2
        window = extract_cleavage_window(protein_seq, 2, 4)
        assert window == "__ABCDEF"
        assert len(window) == 8

    def test_extract_cleavage_window_near_end(self):
        from cleavage_site_profiler import extract_cleavage_window
        protein_seq = "ABCDEFGHIJ"
        # Cleavage at position 8
        window = extract_cleavage_window(protein_seq, 8, 4)
        assert window == "EFGHIJ__"
        assert len(window) == 8

    def test_compute_position_frequencies(self):
        from cleavage_site_profiler import compute_position_frequencies
        windows = ["AABBCCDD", "AABBCCDD"]
        freq = compute_position_frequencies(windows, window=4)
        assert "P4" in freq
        assert "P1" in freq
        assert "P1'" in freq
        assert "P4'" in freq
        assert freq["P4"]["A"] == 1.0

    def test_compute_position_frequencies_empty(self):
        from cleavage_site_profiler import compute_position_frequencies
        assert compute_position_frequencies([], window=4) == {}

    def test_process_neo_nterm_peptides(self):
        from cleavage_site_profiler import load_fasta, process_neo_nterm_peptides

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = self._create_fasta(tmpdir, {"P1": "ABCDEFGHIJKLMNOP"})
            proteins = load_fasta(fasta_path)
            rows = [{"peptide": "GHIJ", "protein": "P1"}]
            result_rows, valid_windows = process_neo_nterm_peptides(rows, proteins, 4)
            assert len(result_rows) == 1
            assert result_rows[0]["cleavage_found"] == "YES"
            assert len(valid_windows) == 1

    def test_process_missing_protein(self):
        from cleavage_site_profiler import process_neo_nterm_peptides
        rows = [{"peptide": "ABCD", "protein": "MISSING"}]
        result_rows, valid_windows = process_neo_nterm_peptides(rows, {}, 4)
        assert result_rows[0]["cleavage_found"] == "NO"
        assert len(valid_windows) == 0

    def test_full_pipeline(self):
        from cleavage_site_profiler import (
            compute_position_frequencies,
            load_fasta,
            process_neo_nterm_peptides,
            write_output,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = self._create_fasta(tmpdir, {"P1": "ABCDEFGHIJKLMNOP"})
            proteins = load_fasta(fasta_path)
            rows = [{"peptide": "GHIJ", "protein": "P1"}]
            result_rows, valid_windows = process_neo_nterm_peptides(rows, proteins, 4)
            freq = compute_position_frequencies(valid_windows, 4)
            output_path = os.path.join(tmpdir, "profile.tsv")
            write_output(output_path, result_rows, freq)
            assert os.path.exists(output_path)
