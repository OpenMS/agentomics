"""Tests for contaminant_database_merger."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


def _create_test_fasta(path, entries_data):
    import pyopenms as oms

    entries = []
    for acc, seq in entries_data:
        e = oms.FASTAEntry()
        e.identifier = acc
        e.sequence = seq
        e.description = ""
        entries.append(e)
    fasta_file = oms.FASTAFile()
    fasta_file.store(path, entries)


def test_get_builtin_contaminants():
    from contaminant_database_merger import get_builtin_contaminants

    contaminants = get_builtin_contaminants("CONT_")
    assert len(contaminants) == 20
    for c in contaminants:
        assert c.identifier.startswith("CONT_")


def test_merge_with_crap():
    from contaminant_database_merger import merge_contaminants

    with tempfile.TemporaryDirectory() as tmp:
        target_path = os.path.join(tmp, "target.fasta")
        output_path = os.path.join(tmp, "merged.fasta")

        _create_test_fasta(target_path, [("P99999", "ACDEFGHIK"), ("Q11111", "MNPQRSTWY")])

        stats = merge_contaminants(target_path, output_path, add_crap=True, prefix="CONT_")
        assert stats["target_count"] == 2
        assert stats["contaminant_count"] == 20
        assert stats["deduplicated_count"] == 22


def test_deduplication():
    import pyopenms as oms
    from contaminant_database_merger import deduplicate

    entries = []
    for acc in ["P12345", "P12345", "P67890"]:
        e = oms.FASTAEntry()
        e.identifier = acc
        e.sequence = "ACDEFGHIK"
        e.description = ""
        entries.append(e)

    result = deduplicate(entries)
    assert len(result) == 2


def test_merge_custom_contaminants():
    from contaminant_database_merger import merge_contaminants

    with tempfile.TemporaryDirectory() as tmp:
        target_path = os.path.join(tmp, "target.fasta")
        contam_path = os.path.join(tmp, "contaminants.fasta")
        output_path = os.path.join(tmp, "merged.fasta")

        _create_test_fasta(target_path, [("P99999", "ACDEFGHIK")])
        _create_test_fasta(contam_path, [("CUSTOM1", "MNPQRST")])

        stats = merge_contaminants(
            target_path, output_path, add_crap=False, contaminants_path=contam_path, prefix="CONT_"
        )
        assert stats["target_count"] == 1
        assert stats["contaminant_count"] == 1
        assert stats["deduplicated_count"] == 2
