"""Tests for hdx_deuterium_uptake."""

import os
import tempfile

from conftest import requires_pyopenms


@requires_pyopenms
class TestHdxDeuteriumUptake:
    def test_count_exchangeable_amides(self):
        from hdx_deuterium_uptake import count_exchangeable_amides
        # PEPTIDEK: 8 residues, 0 prolines, exchangeable = 8 - 0 - 2 = 6
        assert count_exchangeable_amides("PEPTIDEK") == 5  # 1 proline at position 3

    def test_count_exchangeable_amides_with_proline(self):
        from hdx_deuterium_uptake import count_exchangeable_amides
        # PPPPAAAA: 8 residues, 2 prolines at pos>=2, exchangeable = 8 - 2 - 2 = 4
        assert count_exchangeable_amides("PPPPAAAA") == 4

    def test_count_exchangeable_amides_short(self):
        from hdx_deuterium_uptake import count_exchangeable_amides
        # AA: 2 residues, 0 prolines, exchangeable = max(0, 2-0-2) = 0
        assert count_exchangeable_amides("AA") == 0

    def test_get_peptide_mass(self):
        from hdx_deuterium_uptake import get_peptide_mass
        mass = get_peptide_mass("PEPTIDEK")
        assert 927.0 < mass < 928.0

    def test_get_molecular_formula(self):
        from hdx_deuterium_uptake import get_molecular_formula
        formula = get_molecular_formula("PEPTIDEK")
        assert "C" in formula
        assert "H" in formula

    def test_compute_mass_shift(self):
        from hdx_deuterium_uptake import compute_mass_shift
        assert abs(compute_mass_shift(930.0, 928.0) - 2.0) < 1e-6

    def test_compute_fractional_uptake(self):
        from hdx_deuterium_uptake import DEUTERIUM_MASS_SHIFT, compute_fractional_uptake
        # 3 exchangeable, mass shift = 3 * DEUTERIUM_MASS_SHIFT -> 100% uptake
        shift = 3 * DEUTERIUM_MASS_SHIFT
        frac = compute_fractional_uptake(shift, 3)
        assert abs(frac - 1.0) < 1e-6

    def test_compute_fractional_uptake_with_backexchange(self):
        from hdx_deuterium_uptake import DEUTERIUM_MASS_SHIFT, compute_fractional_uptake
        shift = 3 * DEUTERIUM_MASS_SHIFT * 0.8  # 80% of max after 20% back-exchange
        frac = compute_fractional_uptake(shift, 3, back_exchange_fraction=0.2)
        assert abs(frac - 1.0) < 1e-4

    def test_compute_fractional_uptake_zero_exchangeable(self):
        from hdx_deuterium_uptake import compute_fractional_uptake
        assert compute_fractional_uptake(1.0, 0) == 0.0

    def test_compute_uptake_for_peptide(self):
        from hdx_deuterium_uptake import compute_uptake_for_peptide
        result = compute_uptake_for_peptide(
            "AAALAAAK",
            800.0,
            {"10": 802.0, "60": 804.0},
        )
        assert result["sequence"] == "AAALAAAK"
        assert "mass_shift_t10" in result
        assert "fractional_uptake_t60" in result
        assert result["mass_shift_t10"] == 2.0

    def test_group_by_peptide(self):
        from hdx_deuterium_uptake import group_by_peptide
        rows = [
            {"sequence": "PEPTIDEK", "timepoint": "10", "centroid_mass": "930.0"},
            {"sequence": "PEPTIDEK", "timepoint": "60", "centroid_mass": "932.0"},
            {"sequence": "AAALAAAK", "timepoint": "10", "centroid_mass": "702.0"},
        ]
        grouped = group_by_peptide(rows)
        assert len(grouped) == 2
        assert len(grouped["PEPTIDEK"]) == 2

    def test_write_output(self):
        from hdx_deuterium_uptake import write_output
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "uptake.tsv")
            results = [{"sequence": "PEPTIDEK", "mass_shift_t10": 2.0, "fractional_uptake_t10": 0.5}]
            write_output(output_path, results)
            assert os.path.exists(output_path)
