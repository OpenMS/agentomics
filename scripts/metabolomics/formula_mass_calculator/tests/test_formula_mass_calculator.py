"""Tests for formula_mass_calculator."""

from conftest import requires_pyopenms


@requires_pyopenms
class TestFormulaMassCalculator:
    def test_glucose_mass(self):
        from formula_mass_calculator import calculate_formula_mass

        result = calculate_formula_mass("C6H12O6")
        assert abs(result["monoisotopic_mass"] - 180.0634) < 0.001

    def test_neutral_adduct(self):
        from formula_mass_calculator import calculate_formula_mass

        result = calculate_formula_mass("C6H12O6", "[M]")
        assert abs(result["mz"] - result["monoisotopic_mass"]) < 0.001

    def test_mh_adduct(self):
        from formula_mass_calculator import PROTON, calculate_formula_mass

        result = calculate_formula_mass("C6H12O6", "[M+H]+")
        expected_mz = (result["monoisotopic_mass"] + PROTON) / 1
        assert abs(result["mz"] - expected_mz) < 0.001

    def test_sodium_adduct(self):
        from formula_mass_calculator import calculate_formula_mass

        result = calculate_formula_mass("C6H12O6", "[M+Na]+")
        # Na adduct should give higher m/z than H adduct
        result_h = calculate_formula_mass("C6H12O6", "[M+H]+")
        assert result["mz"] > result_h["mz"]

    def test_doubly_charged(self):
        from formula_mass_calculator import calculate_formula_mass

        result = calculate_formula_mass("C6H12O6", "[M+2H]2+")
        assert result["charge"] == 2
        assert result["mz"] < result["monoisotopic_mass"]

    def test_batch_calculate(self):
        from formula_mass_calculator import batch_calculate

        rows = [
            {"formula": "C6H12O6", "adduct": "[M+H]+"},
            {"formula": "C2H6O", "adduct": "[M+Na]+"},
        ]
        results = batch_calculate(rows)
        assert len(results) == 2
        assert results[0]["formula"] == "C6H12O6"
        assert results[1]["formula"] == "C2H6O"

    def test_average_mass(self):
        from formula_mass_calculator import calculate_formula_mass

        result = calculate_formula_mass("C6H12O6")
        assert result["average_mass"] > result["monoisotopic_mass"]
