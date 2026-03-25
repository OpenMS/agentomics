"""Tests for peptide_modification_analyzer."""

from conftest import requires_pyopenms


@requires_pyopenms
class TestPeptideModificationAnalyzer:
    def test_unmodified_peptide(self):
        from peptide_modification_analyzer import analyze_modification

        result = analyze_modification("PEPTIDEK")
        assert result["length"] == 8
        assert len(result["residue_breakdown"]) == 8
        for r in result["residue_breakdown"]:
            assert r["modification"] == ""

    def test_oxidized_methionine(self):
        from peptide_modification_analyzer import analyze_modification

        result = analyze_modification("PEPTM(Oxidation)IDEK")
        residues = result["residue_breakdown"]
        # Position 5 (M) should have Oxidation
        m_residue = residues[4]
        assert m_residue["residue"] == "M"
        assert "Oxidation" in m_residue["modification"]
        assert m_residue["modification_delta_mass"] > 15.0

    def test_mass_consistency(self):
        from peptide_modification_analyzer import analyze_modification

        result = analyze_modification("PEPTIDEK")
        assert result["total_monoisotopic_mass"] > 900

    def test_charge_state(self):
        from peptide_modification_analyzer import PROTON, analyze_modification

        r1 = analyze_modification("PEPTIDEK", charge=1)
        r2 = analyze_modification("PEPTIDEK", charge=2)
        assert r2["mz"] < r1["mz"]
        expected_mz = (r1["total_monoisotopic_mass"] + 2 * PROTON) / 2
        assert abs(r2["mz"] - expected_mz) < 0.001

    def test_residue_positions(self):
        from peptide_modification_analyzer import analyze_modification

        result = analyze_modification("ACDK")
        residues = result["residue_breakdown"]
        assert residues[0]["residue"] == "A"
        assert residues[0]["position"] == 1
        assert residues[3]["residue"] == "K"
        assert residues[3]["position"] == 4
