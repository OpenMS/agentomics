"""Tests for crosslink_mass_calculator."""

import pytest
from conftest import requires_pyopenms


@requires_pyopenms
class TestCrosslinkMassCalculator:
    def test_dss_crosslink(self):
        from crosslink_mass_calculator import crosslinked_mass

        result = crosslinked_mass("PEPTIDEK", "AVLIDR", "DSS", charge=2)
        assert result["crosslinker"] == "DSS"
        assert result["crosslinker_mass"] == pytest.approx(138.068, abs=0.01)
        assert result["total_mass"] > 0
        assert result["charge"] == 2

    def test_dsso_crosslink(self):
        from crosslink_mass_calculator import crosslinked_mass

        result = crosslinked_mass("PEPTIDEK", "AVLIDR", "DSSO", charge=1)
        assert result["crosslinker_mass"] == pytest.approx(158.004, abs=0.01)

    def test_mz_formula(self):
        from crosslink_mass_calculator import PROTON, crosslinked_mass

        result = crosslinked_mass("PEPTIDEK", "AVLIDR", "DSS", charge=3)
        expected_mz = (result["total_mass"] + 3 * PROTON) / 3
        assert result["mz"] == pytest.approx(expected_mz, abs=1e-6)

    def test_custom_crosslinker(self):
        from crosslink_mass_calculator import crosslinked_mass

        result = crosslinked_mass("PEPTIDEK", "AVLIDR", "CUSTOM", charge=1, custom_mass=150.0)
        assert result["crosslinker_mass"] == pytest.approx(150.0)

    def test_unknown_crosslinker_raises(self):
        from crosslink_mass_calculator import crosslinked_mass

        with pytest.raises(ValueError, match="Unknown crosslinker"):
            crosslinked_mass("PEPTIDEK", "AVLIDR", "UNKNOWN", charge=1)

    def test_total_mass_is_sum(self):
        from crosslink_mass_calculator import crosslinked_mass

        result = crosslinked_mass("PEPTIDEK", "AVLIDR", "DSS", charge=1)
        expected = result["mass_peptide1"] + result["mass_peptide2"] + result["crosslinker_mass"]
        assert result["total_mass"] == pytest.approx(expected, abs=1e-6)

    def test_write_tsv(self, tmp_path):
        from crosslink_mass_calculator import crosslinked_mass, write_tsv

        result = crosslinked_mass("PEPTIDEK", "AVLIDR", "DSS", charge=2)
        out = str(tmp_path / "out.tsv")
        write_tsv([result], out)
        with open(out) as fh:
            lines = fh.readlines()
        assert len(lines) == 2  # header + 1 data row
        assert "peptide1" in lines[0]
