"""Tests for inclusion_list_generator."""

import pytest

pytest.importorskip("pyopenms")

from inclusion_list_generator import calculate_mz, generate_inclusion_list


class TestInclusionListGenerator:
    def test_calculate_mz(self):
        mz = calculate_mz("PEPTIDEK", 2)
        assert mz > 0

    def test_charge_affects_mz(self):
        mz1 = calculate_mz("PEPTIDEK", 1)
        mz2 = calculate_mz("PEPTIDEK", 2)
        mz3 = calculate_mz("PEPTIDEK", 3)
        assert mz1 > mz2 > mz3

    def test_thermo_format(self):
        peptides = [{"peptide": "PEPTIDEK"}]
        entries = generate_inclusion_list(peptides, [2, 3], output_format="thermo")
        assert len(entries) == 2
        assert "m/z" in entries[0]
        assert "z" in entries[0]
        assert "Compound" in entries[0]

    def test_generic_format(self):
        peptides = [{"peptide": "PEPTIDEK"}]
        entries = generate_inclusion_list(peptides, [2], output_format="generic")
        assert len(entries) == 1
        assert "mz" in entries[0]
        assert "charge" in entries[0]
        assert "peptide" in entries[0]

    def test_multiple_charges(self):
        peptides = [{"peptide": "PEPTIDEK"}, {"peptide": "TESTPEP"}]
        entries = generate_inclusion_list(peptides, [2, 3, 4], output_format="generic")
        assert len(entries) == 6  # 2 peptides x 3 charges

    def test_with_rt_info(self):
        peptides = [{"peptide": "PEPTIDEK", "rt_start": "10.0", "rt_end": "15.0", "protein": "PROT1"}]
        entries = generate_inclusion_list(peptides, [2], output_format="thermo")
        assert entries[0]["t start (min)"] == "10.0"
        assert entries[0]["t stop (min)"] == "15.0"

    def test_unknown_format(self):
        with pytest.raises(ValueError, match="Unknown format"):
            generate_inclusion_list([{"peptide": "PEP"}], [2], output_format="invalid")

    def test_empty_peptides(self):
        entries = generate_inclusion_list([], [2], output_format="generic")
        assert entries == []
