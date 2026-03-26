"""Tests for decoy_database_generator."""

import pytest

pyopenms = pytest.importorskip("pyopenms")


class TestDecoyDatabaseGenerator:
    def _write_fasta(self, path, entries):
        """Helper to write a FASTA file using pyopenms."""
        import pyopenms as oms

        fasta_entries = []
        for ident, seq in entries:
            entry = oms.FASTAEntry()
            entry.identifier = ident
            entry.description = f"Test protein {ident}"
            entry.sequence = seq
            fasta_entries.append(entry)
        oms.FASTAFile().store(str(path), fasta_entries)

    def test_generate_decoys_reverse(self, tmp_path):
        from decoy_database_generator import generate_decoys

        fasta_in = tmp_path / "target.fasta"
        fasta_out = tmp_path / "target_decoy.fasta"

        proteins = [
            ("sp|P00001|PROT1", "PEPTIDEK"),
            ("sp|P00002|PROT2", "ACDERFGHIK"),
            ("sp|P00003|PROT3", "LMNPQRSTVWY"),
        ]
        self._write_fasta(fasta_in, proteins)

        count = generate_decoys(str(fasta_in), str(fasta_out), method="reverse")
        assert count == 3

        # Verify output has 6 entries (3 target + 3 decoy)
        import pyopenms as oms

        entries = []
        oms.FASTAFile().load(str(fasta_out), entries)
        assert len(entries) == 6

    def test_decoy_prefix_applied(self, tmp_path):
        from decoy_database_generator import generate_decoys

        fasta_in = tmp_path / "target.fasta"
        fasta_out = tmp_path / "target_decoy.fasta"

        proteins = [
            ("sp|P00001|PROT1", "PEPTIDEK"),
            ("sp|P00002|PROT2", "ACDERFGHIK"),
            ("sp|P00003|PROT3", "LMNPQRSTVWY"),
        ]
        self._write_fasta(fasta_in, proteins)

        generate_decoys(str(fasta_in), str(fasta_out), method="reverse")

        import pyopenms as oms

        entries = []
        oms.FASTAFile().load(str(fasta_out), entries)

        target_ids = [e.identifier for e in entries if not e.identifier.startswith("DECOY_")]
        decoy_ids = [e.identifier for e in entries if e.identifier.startswith("DECOY_")]
        assert len(target_ids) == 3
        assert len(decoy_ids) == 3

        # Each target should have a corresponding decoy
        for tid in target_ids:
            assert f"DECOY_{tid}" in decoy_ids

    def test_decoy_sequences_reversed(self, tmp_path):
        from decoy_database_generator import generate_decoys

        fasta_in = tmp_path / "target.fasta"
        fasta_out = tmp_path / "target_decoy.fasta"

        proteins = [
            ("PROT1", "PEPTIDEK"),
        ]
        self._write_fasta(fasta_in, proteins)

        generate_decoys(str(fasta_in), str(fasta_out), method="reverse")

        import pyopenms as oms

        entries = []
        oms.FASTAFile().load(str(fasta_out), entries)

        target_seq = None
        decoy_seq = None
        for e in entries:
            if e.identifier == "PROT1":
                target_seq = e.sequence
            elif e.identifier == "DECOY_PROT1":
                decoy_seq = e.sequence

        assert target_seq is not None
        assert decoy_seq is not None
        # The reversed sequence should be different from the original
        assert target_seq != decoy_seq
        # reverseProtein reverses the entire protein sequence
        assert decoy_seq == target_seq[::-1]

    def test_shuffle_method(self, tmp_path):
        from decoy_database_generator import generate_decoys

        fasta_in = tmp_path / "target.fasta"
        fasta_out = tmp_path / "target_decoy.fasta"

        proteins = [
            ("PROT1", "PEPTIDEKAAAGGGLLLTTTRRR"),
            ("PROT2", "ACDERFGHIKAAAGGGLLLTTTRRR"),
            ("PROT3", "LMNPQRSTVWYAAAGGGLLLTTTRRR"),
        ]
        self._write_fasta(fasta_in, proteins)

        count = generate_decoys(str(fasta_in), str(fasta_out), method="shuffle")
        assert count == 3

        import pyopenms as oms

        entries = []
        oms.FASTAFile().load(str(fasta_out), entries)
        assert len(entries) == 6

        # Decoy sequences should differ from target sequences
        targets = {e.identifier: e.sequence for e in entries if not e.identifier.startswith("DECOY_")}
        decoys = {e.identifier: e.sequence for e in entries if e.identifier.startswith("DECOY_")}
        for tid, tseq in targets.items():
            dseq = decoys.get(f"DECOY_{tid}")
            assert dseq is not None
            # Shuffled sequence should have same length
            assert len(dseq) == len(tseq)

    def test_empty_fasta(self, tmp_path):
        from decoy_database_generator import generate_decoys

        fasta_in = tmp_path / "empty.fasta"
        fasta_out = tmp_path / "empty_decoy.fasta"

        # Create an empty file (pyopenms FASTAFile.store([]) creates empty file)
        fasta_in.write_text("")

        count = generate_decoys(str(fasta_in), str(fasta_out), method="reverse")
        assert count == 0

        # The output should also be an empty file
        assert fasta_out.exists()
        assert fasta_out.stat().st_size == 0
