"""Tests for fasta_subset_extractor."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


def test_filter_by_accessions():
    import pyopenms as oms
    from fasta_subset_extractor import filter_by_accessions

    entries = []
    for acc, seq in [("sp|P12345|PROT1", "ACDEFGHIK"), ("sp|P67890|PROT2", "MNPQRST")]:
        e = oms.FASTAEntry()
        e.identifier = acc
        e.sequence = seq
        e.description = ""
        entries.append(e)

    result = filter_by_accessions(entries, {"P12345"})
    assert len(result) == 1
    assert "P12345" in result[0].identifier


def test_filter_by_keyword():
    import pyopenms as oms
    from fasta_subset_extractor import filter_by_keyword

    e1 = oms.FASTAEntry()
    e1.identifier = "sp|P12345|PROT1"
    e1.description = "Human kinase protein"
    e1.sequence = "ACDEFGHIK"

    e2 = oms.FASTAEntry()
    e2.identifier = "sp|P67890|PROT2"
    e2.description = "Mouse albumin"
    e2.sequence = "MNPQRST"

    result = filter_by_keyword([e1, e2], "kinase")
    assert len(result) == 1


def test_filter_by_length():
    import pyopenms as oms
    from fasta_subset_extractor import filter_by_length

    entries = []
    for seq in ["ACDE", "ACDEFGHIKLMNPQ", "AC"]:
        e = oms.FASTAEntry()
        e.identifier = "test"
        e.sequence = seq
        e.description = ""
        entries.append(e)

    result = filter_by_length(entries, min_length=3, max_length=10)
    assert len(result) == 1
    assert result[0].sequence == "ACDE"


def test_extract_subset_roundtrip():
    import pyopenms as oms
    from fasta_subset_extractor import extract_subset

    # Create a synthetic FASTA file
    entries = []
    for acc, seq in [("sp|P12345|PROT1", "ACDEFGHIK"), ("sp|P67890|PROT2", "MNPQRSTWY")]:
        e = oms.FASTAEntry()
        e.identifier = acc
        e.sequence = seq
        e.description = ""
        entries.append(e)

    with tempfile.TemporaryDirectory() as tmp:
        input_path = os.path.join(tmp, "input.fasta")
        output_path = os.path.join(tmp, "output.fasta")
        accessions_path = os.path.join(tmp, "accessions.txt")

        fasta_file = oms.FASTAFile()
        fasta_file.store(input_path, entries)

        with open(accessions_path, "w") as fh:
            fh.write("P12345\n")

        stats = extract_subset(input_path, output_path, accessions_file=accessions_path)
        assert stats["total_input"] == 2
        assert stats["total_output"] == 1

        result = []
        fasta_file.load(output_path, result)
        assert len(result) == 1
