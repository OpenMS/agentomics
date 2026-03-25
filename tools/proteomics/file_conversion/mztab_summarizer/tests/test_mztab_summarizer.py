"""Tests for mztab_summarizer."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")

SAMPLE_MZTAB = """\
MTD\tmzTab-version\t1.0.0
MTD\tsearch_engine[1]\t[MS, MS:1002995, SEQUEST, ]
PRH\taccession\tdescription\ttaxid\tspecies
PRT\tP12345\tProtein 1\t9606\tHomo sapiens
PRT\tP67890\tProtein 2\t9606\tHomo sapiens
PRT\tP12345\tProtein 1 duplicate\t9606\tHomo sapiens
PEH\tsequence\taccession\tmodifications
PEP\tACDEFGHIK\tP12345\tnull
PEP\tMNPQRSTWY\tP67890\t1-MOD:01214
PEP\tACDEFGHIK\tP12345\tnull
PSH\tsequence\tPSM_ID\taccession
PSM\tACDEFGHIK\t1\tP12345
PSM\tMNPQRSTWY\t2\tP67890
PSM\tACDEFGHIK\t3\tP12345
"""


def test_parse_mztab():
    from mztab_summarizer import parse_mztab

    with tempfile.TemporaryDirectory() as tmp:
        mztab_path = os.path.join(tmp, "test.mzTab")
        with open(mztab_path, "w") as fh:
            fh.write(SAMPLE_MZTAB)

        data = parse_mztab(mztab_path)
        assert len(data["proteins"]) == 3
        assert len(data["peptides"]) == 3
        assert len(data["psms"]) == 3
        assert "mzTab-version" in data["metadata"]


def test_summarize_mztab():
    from mztab_summarizer import summarize_mztab

    with tempfile.TemporaryDirectory() as tmp:
        mztab_path = os.path.join(tmp, "test.mzTab")
        with open(mztab_path, "w") as fh:
            fh.write(SAMPLE_MZTAB)

        summary = summarize_mztab(mztab_path)
        assert summary["protein_count"] == 3
        assert summary["unique_protein_accessions"] == 2
        assert summary["peptide_count"] == 3
        assert summary["unique_peptide_sequences"] == 2
        assert summary["psm_count"] == 3
        assert summary["unique_psm_peptides"] == 2
        assert len(summary["search_engines"]) == 1


def test_write_summary_tsv():
    from mztab_summarizer import summarize_mztab, write_summary_tsv

    with tempfile.TemporaryDirectory() as tmp:
        mztab_path = os.path.join(tmp, "test.mzTab")
        tsv_path = os.path.join(tmp, "summary.tsv")

        with open(mztab_path, "w") as fh:
            fh.write(SAMPLE_MZTAB)

        summary = summarize_mztab(mztab_path)
        write_summary_tsv(summary, tsv_path)

        with open(tsv_path) as fh:
            lines = fh.readlines()
        assert lines[0].strip() == "metric\tvalue"
        assert len(lines) > 1
