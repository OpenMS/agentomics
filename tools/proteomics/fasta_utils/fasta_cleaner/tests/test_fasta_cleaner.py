"""Tests for fasta_cleaner."""

import os
import tempfile

from conftest import requires_pyopenms


def _make_entry(identifier, sequence):
    import pyopenms as oms

    e = oms.FASTAEntry()
    e.identifier = identifier
    e.sequence = sequence
    e.description = ""
    return e


@requires_pyopenms
def test_remove_duplicates():
    from fasta_cleaner import remove_duplicates

    entries = [_make_entry("P1", "ACDEFGHIK"), _make_entry("P2", "ACDEFGHIK"), _make_entry("P3", "MNPQRST")]
    result = remove_duplicates(entries)
    assert len(result) == 2


@requires_pyopenms
def test_remove_stop_codons():
    from fasta_cleaner import remove_stop_codons

    entries = [_make_entry("P1", "ACDEFGHIK*"), _make_entry("P2", "MNPQRST**")]
    result = remove_stop_codons(entries)
    assert result[0].sequence == "ACDEFGHIK"
    assert result[1].sequence == "MNPQRST"


@requires_pyopenms
def test_fix_headers():
    from fasta_cleaner import fix_headers

    entries = [_make_entry("P1  extra   spaces", "ACDE")]
    result = fix_headers(entries)
    assert result[0].identifier == "P1 extra spaces"


@requires_pyopenms
def test_filter_by_length():
    from fasta_cleaner import filter_by_length

    entries = [_make_entry("P1", "AC"), _make_entry("P2", "ACDEFGHIK"), _make_entry("P3", "A" * 100)]
    result = filter_by_length(entries, min_length=3, max_length=50)
    assert len(result) == 1
    assert result[0].identifier == "P2"


@requires_pyopenms
def test_clean_fasta_roundtrip():
    import pyopenms as oms
    from fasta_cleaner import clean_fasta

    with tempfile.TemporaryDirectory() as tmp:
        input_path = os.path.join(tmp, "input.fasta")
        output_path = os.path.join(tmp, "output.fasta")

        entries = [
            _make_entry("P1", "ACDEFGHIK*"),
            _make_entry("P2", "ACDEFGHIK*"),
            _make_entry("P3", "AC"),
        ]
        fasta_file = oms.FASTAFile()
        fasta_file.store(input_path, entries)

        stats = clean_fasta(
            input_path, output_path, dedup=True, min_length=5, strip_stop_codons=True
        )
        assert stats["total_input"] == 3
        assert stats["total_output"] == 1

        result = []
        fasta_file.load(output_path, result)
        assert result[0].sequence == "ACDEFGHIK"
