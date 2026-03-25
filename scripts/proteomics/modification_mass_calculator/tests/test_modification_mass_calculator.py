"""Tests for modification_mass_calculator."""

from conftest import requires_pyopenms


@requires_pyopenms
class TestModificationMassCalculator:
    def test_search_oxidation(self):
        from modification_mass_calculator import search_modification

        results = search_modification("Oxidation")
        assert len(results) > 0
        names = [r["name"] for r in results]
        assert "Oxidation" in names

    def test_search_phospho(self):
        from modification_mass_calculator import search_modification

        results = search_modification("Phospho")
        assert len(results) > 0
        # Phospho should have ~79.966 Da delta mass
        masses = [r["delta_mass"] for r in results]
        assert any(79.9 < m < 80.0 for m in masses)

    def test_list_common_modifications(self):
        from modification_mass_calculator import list_common_modifications

        results = list_common_modifications()
        assert len(results) > 5

    def test_modified_peptide_mass_unmodified(self):
        from modification_mass_calculator import modified_peptide_mass

        result = modified_peptide_mass("PEPTIDEK")
        assert result["monoisotopic_mass"] > 900

    def test_modified_peptide_mass_with_mod(self):
        from modification_mass_calculator import modified_peptide_mass

        unmod = modified_peptide_mass("PEPTMIDEK")
        mod = modified_peptide_mass("PEPTMIDEK", "Oxidation(M):5")
        # Oxidation adds ~16 Da
        assert mod["monoisotopic_mass"] > unmod["monoisotopic_mass"]
        diff = mod["monoisotopic_mass"] - unmod["monoisotopic_mass"]
        assert 15.9 < diff < 16.1

    def test_charge_state(self):
        from modification_mass_calculator import modified_peptide_mass

        r1 = modified_peptide_mass("PEPTIDEK", charge=1)
        r2 = modified_peptide_mass("PEPTIDEK", charge=2)
        assert r2["mz"] < r1["mz"]
