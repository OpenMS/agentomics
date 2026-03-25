"""Tests for fasta_merger."""

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
def test_merge_basic():
    from fasta_merger import merge_fasta_files

    with tempfile.TemporaryDirectory() as tmp:
        f1 = os.path.join(tmp, "db1.fasta")
        f2 = os.path.join(tmp, "db2.fasta")
        out = os.path.join(tmp, "merged.fasta")

        _create_fasta(f1, [("P1", "ACDEFGHIK")])
        _create_fasta(f2, [("P2", "MNPQRSTWY")])

        stats = merge_fasta_files([f1, f2], out)
        assert stats["total_output"] == 2


@requires_pyopenms
def test_merge_dedup_identifier():
    from fasta_merger import merge_fasta_files

    with tempfile.TemporaryDirectory() as tmp:
        f1 = os.path.join(tmp, "db1.fasta")
        f2 = os.path.join(tmp, "db2.fasta")
        out = os.path.join(tmp, "merged.fasta")

        _create_fasta(f1, [("P1", "ACDEFGHIK")])
        _create_fasta(f2, [("P1", "MNPQRSTWY"), ("P2", "ACDEFGHIK")])

        stats = merge_fasta_files([f1, f2], out, remove_duplicates=True, dedup_by="identifier")
        assert stats["total_output"] == 2


@requires_pyopenms
def test_merge_dedup_sequence():
    from fasta_merger import merge_fasta_files

    with tempfile.TemporaryDirectory() as tmp:
        f1 = os.path.join(tmp, "db1.fasta")
        f2 = os.path.join(tmp, "db2.fasta")
        out = os.path.join(tmp, "merged.fasta")

        _create_fasta(f1, [("P1", "ACDEFGHIK")])
        _create_fasta(f2, [("P2", "ACDEFGHIK")])

        stats = merge_fasta_files([f1, f2], out, remove_duplicates=True, dedup_by="sequence")
        assert stats["total_output"] == 1
