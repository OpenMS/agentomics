"""Tests for phospho_motif_analyzer."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestPhosphoMotifAnalyzer:
    def _create_fasta(self, tmpdir, proteins):
        """Helper to create a FASTA file from a dict of accession->sequence."""
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

    def test_extract_window_center(self):
        from phospho_motif_analyzer import extract_window
        seq = "ABCDEFGHIJKLMNOP"
        window = extract_window(seq, 8, 3)  # 1-based pos 8 = 'H'
        assert window == "EFGHIJK"
        assert len(window) == 7

    def test_extract_window_near_start(self):
        from phospho_motif_analyzer import extract_window
        seq = "ABCDEFGHIJ"
        window = extract_window(seq, 2, 3)  # 'B'
        assert window == "__ABCDE"
        assert len(window) == 7

    def test_extract_window_near_end(self):
        from phospho_motif_analyzer import extract_window
        seq = "ABCDEFGHIJ"
        window = extract_window(seq, 9, 3)  # 'I'
        assert window == "FGHIJ__"
        assert len(window) == 7

    def test_validate_peptide(self):
        from phospho_motif_analyzer import validate_peptide
        assert validate_peptide("PEPTIDEK") is True

    def test_extract_motif_windows(self):
        from phospho_motif_analyzer import extract_motif_windows
        proteins = {"P1": "ABCDEFGHIJKLMNOPQRSTUVWXYZ"}
        rows = [{"peptide": "DEFGH", "protein": "P1", "site": "10"}]
        result = extract_motif_windows(rows, proteins, window=3)
        assert len(result) == 1
        assert result[0]["motif_window"] == "GHIJKLM"

    def test_extract_motif_missing_protein(self):
        from phospho_motif_analyzer import extract_motif_windows
        proteins = {}
        rows = [{"peptide": "DEFGH", "protein": "MISSING", "site": "5"}]
        result = extract_motif_windows(rows, proteins, window=3)
        assert result[0]["motif_window"] == "_______"

    def test_compute_position_frequencies(self):
        from phospho_motif_analyzer import compute_position_frequencies
        windows = ["ABA", "ACA", "ADA"]
        freq = compute_position_frequencies(windows, window=1)
        assert freq[-1]["A"] == 1.0
        assert freq[1]["A"] == 1.0
        assert abs(freq[0]["B"] - 1 / 3) < 1e-6

    def test_compute_position_frequencies_empty(self):
        from phospho_motif_analyzer import compute_position_frequencies
        freq = compute_position_frequencies([], window=3)
        assert freq == {}

    def test_format_frequencies(self):
        from phospho_motif_analyzer import format_frequencies
        freq = {0: {"A": 0.5, "B": 0.5}, 1: {"C": 1.0}}
        flat = format_frequencies(freq)
        assert len(flat) == 3

    def test_load_fasta(self):
        from phospho_motif_analyzer import load_fasta
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = self._create_fasta(tmpdir, {"P1": "ACDEFGHIKLMNPQRSTVWY"})
            proteins = load_fasta(fasta_path)
            assert "P1" in proteins
            assert proteins["P1"] == "ACDEFGHIKLMNPQRSTVWY"

    def test_full_pipeline(self):
        from phospho_motif_analyzer import (
            compute_position_frequencies,
            extract_motif_windows,
            load_fasta,
            write_output,
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = self._create_fasta(tmpdir, {"P1": "ACDEFGHIKLMNPQRSTVWY"})
            proteins = load_fasta(fasta_path)
            rows = [{"peptide": "DEFGH", "protein": "P1", "site": "5"}]
            motif_rows = extract_motif_windows(rows, proteins, window=3)
            windows = [r["motif_window"] for r in motif_rows]
            frequencies = compute_position_frequencies(windows, window=3)
            output_path = os.path.join(tmpdir, "output.tsv")
            write_output(output_path, motif_rows, frequencies)
            assert os.path.exists(output_path)
