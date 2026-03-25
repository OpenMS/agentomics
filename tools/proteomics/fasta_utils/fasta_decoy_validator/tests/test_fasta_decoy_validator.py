"""Tests for fasta_decoy_validator."""

import os
import tempfile

from conftest import requires_pyopenms


def _create_fasta(path, entries_data):
    import pyopenms as oms

    entries = []
    for acc, seq in entries_data:
        e = oms.FASTAEntry()
        e.identifier = acc
        e.sequence = seq
        e.description = ""
        entries.append(e)
    oms.FASTAFile().store(path, entries)


@requires_pyopenms
def test_no_decoys():
    from fasta_decoy_validator import validate_decoys

    with tempfile.TemporaryDirectory() as tmp:
        fasta_path = os.path.join(tmp, "target.fasta")
        _create_fasta(fasta_path, [("P12345", "ACDEFGHIK"), ("P67890", "MNPQRSTWY")])

        result = validate_decoys(fasta_path)
        assert result["has_decoys"] is False
        assert result["target_count"] == 2
        assert result["decoy_count"] == 0


@requires_pyopenms
def test_with_decoys():
    from fasta_decoy_validator import validate_decoys

    with tempfile.TemporaryDirectory() as tmp:
        fasta_path = os.path.join(tmp, "td.fasta")
        _create_fasta(fasta_path, [
            ("P12345", "ACDEFGHIK"),
            ("DECOY_P12345", "KIHGFEDCA"),
        ])

        result = validate_decoys(fasta_path)
        assert result["has_decoys"] is True
        assert result["target_count"] == 1
        assert result["decoy_count"] == 1
        assert result["prefix_consistent"] is True
        assert result["decoy_ratio"] == 1.0


@requires_pyopenms
def test_mixed_prefixes():
    from fasta_decoy_validator import validate_decoys

    with tempfile.TemporaryDirectory() as tmp:
        fasta_path = os.path.join(tmp, "mixed.fasta")
        _create_fasta(fasta_path, [
            ("P12345", "ACDEFGHIK"),
            ("DECOY_P12345", "KIHGFEDCA"),
            ("REV_P67890", "YWTSRQPNM"),
        ])

        result = validate_decoys(fasta_path, decoy_prefix="DECOY_")
        assert result["has_decoys"] is True
        assert result["prefix_consistent"] is False
        assert "REV_" in result["alternative_prefixes"]


@requires_pyopenms
def test_reversed_match():
    from fasta_decoy_validator import validate_decoys

    with tempfile.TemporaryDirectory() as tmp:
        fasta_path = os.path.join(tmp, "rev.fasta")
        seq = "ACDEFGHIK"
        _create_fasta(fasta_path, [
            ("P12345", seq),
            ("DECOY_P12345", seq[::-1]),
        ])

        result = validate_decoys(fasta_path)
        assert result["reversed_match_count"] == 1
