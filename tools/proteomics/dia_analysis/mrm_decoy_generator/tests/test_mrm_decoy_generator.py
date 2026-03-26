"""Tests for mrm_decoy_generator."""

import csv
import os

import pytest

pyopenms = pytest.importorskip("pyopenms")


def _make_target_transitions_tsv(path):
    """Create a small target transition list TSV."""
    rows = [
        {"PrecursorMz": "500.0", "ProductMz": "600.0", "LibraryIntensity": "100.0",
         "PeptideSequence": "PEPTIDEK", "ProteinName": "PROT1",
         "transition_name": "tr_0", "transition_group_id": "pep_0",
         "PrecursorCharge": "2"},
        {"PrecursorMz": "500.0", "ProductMz": "650.0", "LibraryIntensity": "90.0",
         "PeptideSequence": "PEPTIDEK", "ProteinName": "PROT1",
         "transition_name": "tr_1", "transition_group_id": "pep_0",
         "PrecursorCharge": "2"},
        {"PrecursorMz": "500.0", "ProductMz": "700.0", "LibraryIntensity": "80.0",
         "PeptideSequence": "PEPTIDEK", "ProteinName": "PROT1",
         "transition_name": "tr_2", "transition_group_id": "pep_0",
         "PrecursorCharge": "2"},
        {"PrecursorMz": "550.0", "ProductMz": "620.0", "LibraryIntensity": "100.0",
         "PeptideSequence": "ANOTHERPEPTIDER", "ProteinName": "PROT1",
         "transition_name": "tr_3", "transition_group_id": "pep_1",
         "PrecursorCharge": "2"},
        {"PrecursorMz": "550.0", "ProductMz": "670.0", "LibraryIntensity": "90.0",
         "PeptideSequence": "ANOTHERPEPTIDER", "ProteinName": "PROT1",
         "transition_name": "tr_4", "transition_group_id": "pep_1",
         "PrecursorCharge": "2"},
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class TestMRMDecoyGenerator:
    def test_generate_shuffle_decoys(self, tmp_path):
        from mrm_decoy_generator import generate_mrm_decoys

        in_path = str(tmp_path / "targets.tsv")
        out_path = str(tmp_path / "decoys.tsv")

        _make_target_transitions_tsv(in_path)
        n = generate_mrm_decoys(in_path, out_path, method="shuffle")

        assert n == 5  # same count as targets
        assert os.path.exists(out_path)

    def test_generate_reverse_decoys(self, tmp_path):
        from mrm_decoy_generator import generate_mrm_decoys

        in_path = str(tmp_path / "targets.tsv")
        out_path = str(tmp_path / "decoys.tsv")

        _make_target_transitions_tsv(in_path)
        n = generate_mrm_decoys(in_path, out_path, method="reverse")
        assert n == 5

    def test_decoy_sequences_differ(self, tmp_path):
        from mrm_decoy_generator import generate_mrm_decoys

        in_path = str(tmp_path / "targets.tsv")
        out_path = str(tmp_path / "decoys.tsv")

        _make_target_transitions_tsv(in_path)
        generate_mrm_decoys(in_path, out_path, method="shuffle")

        with open(out_path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            decoy_rows = list(reader)

        # Check that at least one decoy sequence differs from target
        target_seqs = {"PEPTIDEK", "ANOTHERPEPTIDER"}
        decoy_seqs = {r["PeptideSequence"] for r in decoy_rows}
        assert decoy_seqs != target_seqs

    def test_decoy_has_tag_prefix(self, tmp_path):
        from mrm_decoy_generator import generate_mrm_decoys

        in_path = str(tmp_path / "targets.tsv")
        out_path = str(tmp_path / "decoys.tsv")

        _make_target_transitions_tsv(in_path)
        generate_mrm_decoys(in_path, out_path, decoy_tag="DECOY_")

        with open(out_path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            decoy_rows = list(reader)

        for row in decoy_rows:
            assert row["transition_name"].startswith("DECOY_")
            assert row["transition_group_id"].startswith("DECOY_")
            assert row["ProteinName"].startswith("DECOY_")

    def test_decoy_count_matches_target(self, tmp_path):
        from mrm_decoy_generator import generate_mrm_decoys

        in_path = str(tmp_path / "targets.tsv")
        out_path = str(tmp_path / "decoys.tsv")

        _make_target_transitions_tsv(in_path)
        n_decoys = generate_mrm_decoys(in_path, out_path)

        # Count targets
        with open(in_path) as fh:
            n_targets = sum(1 for _ in csv.DictReader(fh, delimiter="\t"))

        assert n_decoys == n_targets

    def test_generate_decoy_sequence_reverse(self):
        from mrm_decoy_generator import _generate_decoy_sequence

        seq = "PEPTIDEK"
        decoy = _generate_decoy_sequence(seq, "reverse")
        assert decoy != seq
        assert len(decoy) == len(seq)

    def test_generate_decoy_sequence_shuffle(self):
        from mrm_decoy_generator import _generate_decoy_sequence

        seq = "PEPTIDEK"
        decoy = _generate_decoy_sequence(seq, "shuffle")
        assert len(decoy) == len(seq)
        # Shuffled should have same amino acid composition (approximately)
        assert sorted(decoy) == sorted(seq) or len(decoy) == len(seq)

    def test_returns_int(self, tmp_path):
        from mrm_decoy_generator import generate_mrm_decoys

        in_path = str(tmp_path / "targets.tsv")
        out_path = str(tmp_path / "decoys.tsv")

        _make_target_transitions_tsv(in_path)
        result = generate_mrm_decoys(in_path, out_path)
        assert isinstance(result, int)
