"""Tests for phospho_enrichment_qc."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestPhosphoEnrichmentQC:
    def test_is_phosphopeptide_true(self):
        from phospho_enrichment_qc import is_phosphopeptide
        assert is_phosphopeptide("PEPTIDES(Phospho)K") is True
        assert is_phosphopeptide("PEPTIDET[80]K") is True

    def test_is_phosphopeptide_false(self):
        from phospho_enrichment_qc import is_phosphopeptide
        assert is_phosphopeptide("PEPTIDEK") is False

    def test_count_phospho_residues_ser(self):
        from phospho_enrichment_qc import count_phospho_residues
        counts = count_phospho_residues("PEPTIDES(Phospho)K")
        assert counts["pSer"] == 1
        assert counts["pThr"] == 0
        assert counts["pTyr"] == 0

    def test_count_phospho_residues_thr(self):
        from phospho_enrichment_qc import count_phospho_residues
        counts = count_phospho_residues("PEPTIDET(Phospho)K")
        assert counts["pThr"] == 1

    def test_count_phospho_residues_tyr(self):
        from phospho_enrichment_qc import count_phospho_residues
        counts = count_phospho_residues("PEPTIDEY(Phospho)K")
        assert counts["pTyr"] == 1

    def test_count_phospho_multiple(self):
        from phospho_enrichment_qc import count_phospho_residues
        counts = count_phospho_residues("S(Phospho)EPTIDES(Phospho)T(Phospho)K")
        assert counts["pSer"] == 2
        assert counts["pThr"] == 1

    def test_get_peptide_length(self):
        from phospho_enrichment_qc import get_peptide_length
        assert get_peptide_length("PEPTIDEK") == 8

    def test_compute_enrichment_stats(self):
        from phospho_enrichment_qc import compute_enrichment_stats
        rows = [
            {"sequence": "PEPTIDES(Phospho)K"},
            {"sequence": "PEPTIDET(Phospho)K"},
            {"sequence": "PEPTIDEK"},
            {"sequence": "ACDEFGHIK"},
        ]
        counts, ratios = compute_enrichment_stats(rows)
        assert counts["total"] == 4
        assert counts["phospho"] == 2
        assert counts["pSer"] == 1
        assert counts["pThr"] == 1
        assert abs(ratios["enrichment_efficiency"] - 0.5) < 1e-6

    def test_compute_enrichment_stats_empty(self):
        from phospho_enrichment_qc import compute_enrichment_stats
        counts, ratios = compute_enrichment_stats([])
        assert counts["total"] == 0
        assert ratios["enrichment_efficiency"] == 0.0

    def test_write_output(self):
        from phospho_enrichment_qc import write_output
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "enrichment.tsv")
            counts = {"total": 10, "phospho": 7, "pSer": 5, "pThr": 1, "pTyr": 1}
            ratios = {
                "enrichment_efficiency": 0.7,
                "pSer_ratio": 5 / 7, "pThr_ratio": 1 / 7, "pTyr_ratio": 1 / 7
            }
            write_output(output_path, counts, ratios)
            assert os.path.exists(output_path)
            with open(output_path) as f:
                content = f.read()
            assert "enrichment_efficiency" in content
