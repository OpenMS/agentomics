"""Tests for missed_cleavage_analyzer."""

import pytest

pytest.importorskip("pyopenms")


class TestMissedCleavageAnalyzer:
    def test_no_missed_cleavages(self):
        from missed_cleavage_analyzer import count_missed_cleavages

        # Trypsin cleaves after K/R; a peptide ending in K with no internal K/R
        mc = count_missed_cleavages("PEPTIDEK", "Trypsin")
        assert mc == 0

    def test_one_missed_cleavage(self):
        from missed_cleavage_analyzer import count_missed_cleavages

        # Internal K followed by more residues then ending in K
        mc = count_missed_cleavages("PEPKIDEK", "Trypsin")
        assert mc == 1

    def test_two_missed_cleavages(self):
        from missed_cleavage_analyzer import count_missed_cleavages

        mc = count_missed_cleavages("PEPKIDRKEK", "Trypsin")
        assert mc >= 2

    def test_analyze_distribution(self):
        from missed_cleavage_analyzer import analyze_missed_cleavages

        peptides = ["PEPTIDEK", "PEPKIDEK", "ANOTHERSEQ"]
        analysis = analyze_missed_cleavages(peptides, "Trypsin")
        assert analysis["total_peptides"] == 3
        assert analysis["enzyme"] == "Trypsin"
        assert "distribution" in analysis
        assert len(analysis["peptide_results"]) == 3

    def test_average_mc(self):
        from missed_cleavage_analyzer import analyze_missed_cleavages

        peptides = ["PEPTIDEK", "PEPTIDEK"]
        analysis = analyze_missed_cleavages(peptides, "Trypsin")
        assert analysis["average_missed_cleavages"] == 0.0

    def test_empty_list(self):
        from missed_cleavage_analyzer import analyze_missed_cleavages

        analysis = analyze_missed_cleavages([], "Trypsin")
        assert analysis["total_peptides"] == 0
