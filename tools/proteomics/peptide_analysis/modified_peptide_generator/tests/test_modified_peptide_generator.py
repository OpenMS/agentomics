"""Tests for modified_peptide_generator."""

import pytest

pytest.importorskip("pyopenms")


class TestModifiedPeptideGenerator:
    def test_no_mods(self):
        from modified_peptide_generator import generate_variants

        variants = generate_variants("PEPTIDEK", [], [])
        assert len(variants) == 1
        assert variants[0]["modifications"] == "none"

    def test_oxidation_on_methionine(self):
        from modified_peptide_generator import generate_variants

        variants = generate_variants("PEPTMIDEK", ["Oxidation"], max_mods=1)
        # Should have unmodified + 1 oxidized variant
        assert len(variants) == 2
        masses = [v["monoisotopic_mass"] for v in variants]
        assert masses[1] > masses[0]  # Oxidation adds mass

    def test_max_mods_limit(self):
        from modified_peptide_generator import generate_variants

        # Two M residues, max 1 mod
        variants = generate_variants("MMPEPTIDEK", ["Oxidation"], max_mods=1)
        max_mod_count = max(v["num_modifications"] for v in variants)
        assert max_mod_count <= 1

    def test_find_modifiable_sites(self):
        from modified_peptide_generator import find_modifiable_sites

        sites = find_modifiable_sites("PEPTMIDEK", "Oxidation")
        assert len(sites) == 1
        assert sites[0] == (5, "M")

    def test_fixed_mods(self):
        from modified_peptide_generator import generate_variants

        variants = generate_variants("CPEPTIDEK", [], fixed_mods=["Carbamidomethyl"])
        assert len(variants) == 1
        assert variants[0]["num_modifications"] == 1

    def test_multiple_variable_mods(self):
        from modified_peptide_generator import generate_variants

        variants = generate_variants("MSPEPTIDEK", ["Oxidation", "Phospho"], max_mods=2)
        # M can be oxidized, S can be phosphorylated
        assert len(variants) >= 3  # unmodified, ox-only, phos-only, both
