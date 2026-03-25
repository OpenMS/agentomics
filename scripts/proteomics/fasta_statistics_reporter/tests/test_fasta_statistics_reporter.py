"""Tests for fasta_statistics_reporter."""

import json
import os
import tempfile

from conftest import requires_pyopenms


def _create_test_fasta(path):
    import pyopenms as oms

    entries = []
    for acc, seq in [
        ("sp|P12345|PROT1", "ACDEFGHIKLMNPQR"),
        ("sp|P67890|PROT2", "ACDEFGHIK"),
        ("sp|P12345|PROT1", "ACDEFGHIKLMNPQR"),  # duplicate
    ]:
        e = oms.FASTAEntry()
        e.identifier = acc
        e.sequence = seq
        e.description = ""
        entries.append(e)
    fasta_file = oms.FASTAFile()
    fasta_file.store(path, entries)


@requires_pyopenms
def test_compute_statistics_basic():
    from fasta_statistics_reporter import compute_statistics

    with tempfile.TemporaryDirectory() as tmp:
        fasta_path = os.path.join(tmp, "test.fasta")
        _create_test_fasta(fasta_path)

        stats = compute_statistics(fasta_path)
        assert stats["protein_count"] == 3
        assert stats["length_stats"]["min"] == 9
        assert stats["length_stats"]["max"] == 15
        assert len(stats["duplicate_accessions"]) == 1


@requires_pyopenms
def test_compute_statistics_with_enzyme():
    from fasta_statistics_reporter import compute_statistics

    with tempfile.TemporaryDirectory() as tmp:
        fasta_path = os.path.join(tmp, "test.fasta")
        _create_test_fasta(fasta_path)

        stats = compute_statistics(fasta_path, enzyme="Trypsin")
        assert "tryptic_peptide_count" in stats
        assert stats["tryptic_peptide_count"] > 0


@requires_pyopenms
def test_aa_frequency():
    import pyopenms as oms
    from fasta_statistics_reporter import compute_aa_frequency

    e = oms.FASTAEntry()
    e.identifier = "test"
    e.sequence = "AAACCC"
    e.description = ""

    freq = compute_aa_frequency([e])
    assert freq["A"] == 3
    assert freq["C"] == 3


@requires_pyopenms
def test_output_json():
    from fasta_statistics_reporter import compute_statistics

    with tempfile.TemporaryDirectory() as tmp:
        fasta_path = os.path.join(tmp, "test.fasta")
        _create_test_fasta(fasta_path)

        stats = compute_statistics(fasta_path)
        output = json.dumps(stats, indent=2)
        parsed = json.loads(output)
        assert parsed["protein_count"] == 3
