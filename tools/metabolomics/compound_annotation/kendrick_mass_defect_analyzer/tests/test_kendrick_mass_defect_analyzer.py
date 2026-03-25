"""Tests for kendrick_mass_defect_analyzer."""

import pytest

pytest.importorskip("pyopenms")


class TestComputeKMD:
    def test_ch2_base(self):
        from kendrick_mass_defect_analyzer import compute_kmd_from_formula

        result = compute_kmd_from_formula("C16H32O2", "CH2")
        assert "kendrick_mass" in result
        assert "kmd" in result
        assert isinstance(result["kmd"], float)

    def test_kmd_homologous_series_similar(self):
        """Fatty acids differing by CH2 should have similar KMD."""
        from kendrick_mass_defect_analyzer import compute_kmd_from_formula

        r1 = compute_kmd_from_formula("C16H32O2", "CH2")
        r2 = compute_kmd_from_formula("C18H36O2", "CH2")
        assert abs(r1["kmd"] - r2["kmd"]) < 0.01

    def test_compute_kmd_by_mass(self):
        from kendrick_mass_defect_analyzer import compute_kmd

        result = compute_kmd(256.2402, "CH2")
        assert "kendrick_mass" in result
        assert "kmd" in result


class TestGetBaseExactMass:
    def test_ch2_mass(self):
        from kendrick_mass_defect_analyzer import get_base_exact_mass

        mass = get_base_exact_mass("CH2")
        assert mass == pytest.approx(14.01565, abs=0.001)

    def test_cf2_mass(self):
        from kendrick_mass_defect_analyzer import get_base_exact_mass

        mass = get_base_exact_mass("CF2")
        assert mass > 49.0


class TestGroupHomologousSeries:
    def test_grouping(self):
        from kendrick_mass_defect_analyzer import compute_kmd_from_formula, group_homologous_series

        formulas = ["C14H28O2", "C16H32O2", "C18H36O2", "C6H12O6"]
        results = [compute_kmd_from_formula(f, "CH2") for f in formulas]
        grouped = group_homologous_series(results, kmd_tolerance=0.01)
        assert all("series" in r for r in grouped)
        # Fatty acids should be in same series
        fatty_acid_series = set()
        for r in grouped:
            if r["formula"] in ("C14H28O2", "C16H32O2", "C18H36O2"):
                fatty_acid_series.add(r["series"])
        assert len(fatty_acid_series) == 1

    def test_empty_input(self):
        from kendrick_mass_defect_analyzer import group_homologous_series

        assert group_homologous_series([]) == []
