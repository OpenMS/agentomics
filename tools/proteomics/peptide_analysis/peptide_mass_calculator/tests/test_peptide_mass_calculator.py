"""Tests for peptide_mass_calculator."""


from conftest import requires_pyopenms


@requires_pyopenms
class TestPeptideMassCalculator:
    def test_basic_mass(self):
        from peptide_mass_calculator import peptide_masses

        result = peptide_masses("PEPTIDEK")
        assert result["sequence"] == "PEPTIDEK"
        assert result["charge"] == 1
        assert 927.0 < result["monoisotopic_mass"] < 928.0
        assert result["mz_monoisotopic"] > result["monoisotopic_mass"]

    def test_charge_state(self):
        from peptide_mass_calculator import peptide_masses

        r1 = peptide_masses("PEPTIDEK", charge=1)
        r2 = peptide_masses("PEPTIDEK", charge=2)
        assert r2["mz_monoisotopic"] < r1["mz_monoisotopic"]

    def test_fragment_ions(self):
        from peptide_mass_calculator import fragment_ions

        ions = fragment_ions("PEPTIDEK")
        seq_len = len("PEPTIDEK")
        assert len(ions["b_ions"]) == seq_len - 1
        assert len(ions["y_ions"]) == seq_len - 1

    def test_modified_sequence(self):
        from peptide_mass_calculator import peptide_masses

        result = peptide_masses("PEPTM[147]IDEK")
        assert result["monoisotopic_mass"] > 0

    def test_mz_formula(self):
        from peptide_mass_calculator import PROTON, peptide_masses

        r = peptide_masses("PEPTIDEK", charge=2)
        expected = (r["monoisotopic_mass"] + 2 * PROTON) / 2
        assert abs(r["mz_monoisotopic"] - expected) < 1e-6
