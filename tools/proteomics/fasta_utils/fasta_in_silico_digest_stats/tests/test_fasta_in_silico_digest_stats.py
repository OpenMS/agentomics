"""Tests for fasta_in_silico_digest_stats."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


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


def test_digest_basic():
    from fasta_in_silico_digest_stats import digest_fasta

    with tempfile.TemporaryDirectory() as tmp:
        fasta_path = os.path.join(tmp, "test.fasta")
        # Two tryptic cleavage sites (after K and R)
        _create_fasta(fasta_path, [("P1", "ACDEFGHIKLMNPQRSTWYACDEFGHIK")])

        stats = digest_fasta(fasta_path, enzyme="Trypsin", missed_cleavages=0, min_length=6)
        assert stats["protein_count"] == 1
        assert stats["total_peptides"] > 0
        assert stats["unique_peptides"] > 0


def test_digest_with_missed_cleavages():
    from fasta_in_silico_digest_stats import digest_fasta

    with tempfile.TemporaryDirectory() as tmp:
        fasta_path = os.path.join(tmp, "test.fasta")
        _create_fasta(fasta_path, [("P1", "ACDEFGHIKLMNPQRSTWYACDEFGHIK")])

        stats_0 = digest_fasta(fasta_path, enzyme="Trypsin", missed_cleavages=0, min_length=6)
        stats_2 = digest_fasta(fasta_path, enzyme="Trypsin", missed_cleavages=2, min_length=6)
        assert stats_2["total_peptides"] >= stats_0["total_peptides"]


def test_write_tsv():
    from fasta_in_silico_digest_stats import digest_fasta, write_tsv

    with tempfile.TemporaryDirectory() as tmp:
        fasta_path = os.path.join(tmp, "test.fasta")
        tsv_path = os.path.join(tmp, "stats.tsv")
        _create_fasta(fasta_path, [("P1", "ACDEFGHIKLMNPQRSTWYACDEFGHIK")])

        stats = digest_fasta(fasta_path, enzyme="Trypsin", min_length=6)
        write_tsv(stats, tsv_path)

        with open(tsv_path) as fh:
            lines = fh.readlines()
        assert lines[0].strip().startswith("sequence")
        assert len(lines) > 1


def test_mass_stats():
    from fasta_in_silico_digest_stats import digest_fasta

    with tempfile.TemporaryDirectory() as tmp:
        fasta_path = os.path.join(tmp, "test.fasta")
        _create_fasta(fasta_path, [("P1", "ACDEFGHIKLMNPQRSTWYACDEFGHIK")])

        stats = digest_fasta(fasta_path, enzyme="Trypsin", min_length=6)
        assert "min" in stats["mass_stats"]
        assert "max" in stats["mass_stats"]
        assert stats["mass_stats"]["min"] > 0
