"""Tests for fasta_taxonomy_splitter."""

import os
import tempfile

from conftest import requires_pyopenms


def _create_fasta(path, entries_data):
    """Create FASTA with (identifier, description, sequence) tuples."""
    import pyopenms as oms

    entries = []
    for identifier, desc, seq in entries_data:
        e = oms.FASTAEntry()
        e.identifier = identifier
        e.description = desc
        e.sequence = seq
        entries.append(e)
    oms.FASTAFile().store(path, entries)


@requires_pyopenms
def test_split_by_taxonomy():
    from fasta_taxonomy_splitter import split_by_taxonomy

    with tempfile.TemporaryDirectory() as tmp:
        fasta_path = os.path.join(tmp, "combined.fasta")
        out_dir = os.path.join(tmp, "split")

        _create_fasta(fasta_path, [
            ("sp|P12345|PROT1", "Protein1 OS=Homo sapiens OX=9606", "ACDEFGHIK"),
            ("sp|P67890|PROT2", "Protein2 OS=Homo sapiens OX=9606", "MNPQRSTWY"),
            ("sp|Q11111|PROT3", "Protein3 OS=Mus musculus OX=10090", "ACDEFGHIK"),
        ])

        stats = split_by_taxonomy(fasta_path, out_dir, pattern=r"OS=([^=]+)\s+OX=")
        assert stats["total_entries"] == 3
        assert stats["taxonomy_groups"] == 2
        assert "Homo sapiens" in stats["files"]
        assert stats["files"]["Homo sapiens"]["count"] == 2


@requires_pyopenms
def test_unmatched_entries():
    from fasta_taxonomy_splitter import split_by_taxonomy

    with tempfile.TemporaryDirectory() as tmp:
        fasta_path = os.path.join(tmp, "combined.fasta")
        out_dir = os.path.join(tmp, "split")

        _create_fasta(fasta_path, [
            ("sp|P12345|PROT1", "Protein1 OS=Homo sapiens OX=9606", "ACDEFGHIK"),
            ("CUSTOM_PROT", "No taxonomy info", "MNPQRSTWY"),
        ])

        stats = split_by_taxonomy(fasta_path, out_dir, pattern=r"OS=([^=]+)\s+OX=")
        assert stats["unmatched_count"] == 1


@requires_pyopenms
def test_sanitize_filename():
    from fasta_taxonomy_splitter import sanitize_filename

    assert sanitize_filename("Homo sapiens") == "Homo_sapiens"
    assert sanitize_filename("E. coli (strain K12)") == "E_coli_strain_K12"
