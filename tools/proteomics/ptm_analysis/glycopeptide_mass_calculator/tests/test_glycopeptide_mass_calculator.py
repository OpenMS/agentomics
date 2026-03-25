"""Tests for glycopeptide_mass_calculator."""

import pytest

pytest.importorskip("pyopenms")


class TestGlycopeptideMassCalculator:
    def test_basic_glycopeptide(self):
        from glycopeptide_mass_calculator import glycopeptide_mass

        result = glycopeptide_mass("PEPTIDEK", "HexNAc(2)Hex(5)", charge=1)
        assert result["sequence"] == "PEPTIDEK"
        assert result["glycan_mass"] > 0
        assert result["total_mass"] > result["peptide_mass"]

    def test_glycan_mass_calculation(self):
        from glycopeptide_mass_calculator import glycan_mass, parse_glycan

        comp = parse_glycan("HexNAc(2)Hex(5)Fuc(1)")
        mass = glycan_mass(comp)
        expected = 203.079 * 2 + 162.053 * 5 + 146.058 * 1
        assert mass == pytest.approx(expected, abs=0.01)

    def test_parse_glycan(self):
        from glycopeptide_mass_calculator import parse_glycan

        comp = parse_glycan("HexNAc(2)Hex(3)NeuAc(1)")
        assert comp == {"HexNAc": 2, "Hex": 3, "NeuAc": 1}

    def test_parse_glycan_invalid(self):
        from glycopeptide_mass_calculator import parse_glycan

        with pytest.raises(ValueError, match="Could not parse"):
            parse_glycan("invalid")

    def test_parse_glycan_unknown_residue(self):
        from glycopeptide_mass_calculator import parse_glycan

        with pytest.raises(ValueError, match="Unknown glycan residue"):
            parse_glycan("Unknown(3)")

    def test_mz_formula(self):
        from glycopeptide_mass_calculator import PROTON, glycopeptide_mass

        result = glycopeptide_mass("PEPTIDEK", "HexNAc(2)Hex(5)", charge=3)
        expected_mz = (result["total_mass"] + 3 * PROTON) / 3
        assert result["mz"] == pytest.approx(expected_mz, abs=1e-6)

    def test_total_mass_is_sum(self):
        from glycopeptide_mass_calculator import glycopeptide_mass

        result = glycopeptide_mass("PEPTIDEK", "HexNAc(2)Hex(5)", charge=1)
        expected = result["peptide_mass"] + result["glycan_mass"]
        assert result["total_mass"] == pytest.approx(expected, abs=1e-6)

    def test_write_tsv(self, tmp_path):
        from glycopeptide_mass_calculator import glycopeptide_mass, write_tsv

        result = glycopeptide_mass("PEPTIDEK", "HexNAc(2)Hex(5)", charge=2)
        out = str(tmp_path / "out.tsv")
        write_tsv([result], out)
        with open(out) as fh:
            lines = fh.readlines()
        assert len(lines) == 2
        assert "sequence" in lines[0]
