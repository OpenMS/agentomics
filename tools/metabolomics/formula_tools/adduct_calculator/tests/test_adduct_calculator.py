"""Tests for adduct_calculator."""

import pytest

pytest.importorskip("pyopenms")


class TestAdductCalculator:
    def test_formula_to_mass(self):
        from adduct_calculator import formula_to_mass

        mass = formula_to_mass("C6H12O6")
        assert abs(mass - 180.0634) < 0.001

    def test_positive_adducts(self):
        from adduct_calculator import compute_adducts

        results = compute_adducts(180.0634, mode="positive")
        assert len(results) > 0
        names = [r["adduct"] for r in results]
        assert "[M+H]+" in names
        assert "[M+Na]+" in names
        assert "[M+K]+" in names

    def test_negative_adducts(self):
        from adduct_calculator import compute_adducts

        results = compute_adducts(180.0634, mode="negative")
        assert len(results) > 0
        names = [r["adduct"] for r in results]
        assert "[M-H]-" in names
        assert "[M+Cl]-" in names

    def test_mh_plus(self):
        from adduct_calculator import PROTON, compute_adducts

        results = compute_adducts(180.0634, mode="positive")
        mh = next(r for r in results if r["adduct"] == "[M+H]+")
        expected = (180.0634 + PROTON) / 1
        assert abs(mh["mz"] - expected) < 0.001

    def test_doubly_charged(self):
        from adduct_calculator import compute_adducts

        results = compute_adducts(180.0634, mode="positive")
        m2h = next(r for r in results if r["adduct"] == "[M+2H]2+")
        assert m2h["charge"] == 2
        assert m2h["mz"] < 180.0634  # Should be about half

    def test_all_mz_positive(self):
        from adduct_calculator import compute_adducts

        for r in compute_adducts(500.0, mode="positive"):
            assert r["mz"] > 0
