"""Tests for protein_group_reporter."""

import pytest

pytest.importorskip("pyopenms")


class TestProteinGroupReporter:
    def test_map_peptides_to_proteins(self):
        from protein_group_reporter import map_peptides_to_proteins

        proteins = {"P1": "PEPTIDEKAVLIDR", "P2": "ACDEFGHIK"}
        peptides = ["PEPTIDEK", "AVLIDR", "ACDEFG"]
        result = map_peptides_to_proteins(peptides, proteins)
        assert "P1" in result
        assert "PEPTIDEK" in result["P1"]
        assert "AVLIDR" in result["P1"]
        assert "P2" in result
        assert "ACDEFG" in result["P2"]

    def test_compute_sequence_coverage(self):
        from protein_group_reporter import compute_sequence_coverage

        protein = "PEPTIDEKAVLIDR"  # 14 aa
        coverage = compute_sequence_coverage(protein, {"PEPTIDEK"})  # 8 aa covered
        assert 0.5 < coverage < 0.6  # 8/14 = 0.571

    def test_compute_coverage_full(self):
        from protein_group_reporter import compute_sequence_coverage

        protein = "PEPTIDEK"
        coverage = compute_sequence_coverage(protein, {"PEPTIDEK"})
        assert coverage == 1.0

    def test_compute_coverage_empty(self):
        from protein_group_reporter import compute_sequence_coverage

        assert compute_sequence_coverage("PEPTIDEK", set()) == 0.0
        assert compute_sequence_coverage("", {"PEPTIDEK"}) == 0.0

    def test_build_protein_groups(self):
        from protein_group_reporter import build_protein_groups

        protein_peptides = {
            "P1": {"PEPTIDEK", "AVLIDR"},
            "P2": {"PEPTIDEK", "AVLIDR"},  # Same set -> same group
            "P3": {"ACDEFG"},
        }
        proteins = {"P1": "PEPTIDEKAVLIDR", "P2": "PEPTIDEKAVLIDR", "P3": "ACDEFGHIK"}
        groups = build_protein_groups(protein_peptides, proteins)
        assert len(groups) == 2
        # The group with more peptides should be first
        assert groups[0]["n_peptides"] == 2

    def test_write_tsv(self, tmp_path):
        from protein_group_reporter import write_tsv

        results = [{
            "protein_group": "P1;P2", "lead_protein": "P1",
            "n_proteins": 2, "n_peptides": 2,
            "peptides": "AVLIDR;PEPTIDEK", "sequence_coverage": 0.571,
        }]
        out = str(tmp_path / "groups.tsv")
        write_tsv(results, out)
        with open(out) as fh:
            lines = fh.readlines()
        assert len(lines) == 2
        assert "protein_group" in lines[0]
