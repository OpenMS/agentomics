"""Tests for mass_accuracy_calculator."""


from conftest import requires_pyopenms


@requires_pyopenms
class TestMassAccuracyCalculator:
    def test_sequence_theoretical(self):
        from mass_accuracy_calculator import theoretical_mz_from_sequence

        mz = theoretical_mz_from_sequence("PEPTIDEK", 1)
        assert 928.0 < mz < 929.0

    def test_formula_theoretical(self):
        from mass_accuracy_calculator import theoretical_mz_from_formula

        mz = theoretical_mz_from_formula("C6H12O6", 1)
        assert 181.0 < mz < 182.0

    def test_ppm_zero_error(self):
        from mass_accuracy_calculator import ppm_error

        assert ppm_error(500.0, 500.0) == 0.0

    def test_ppm_positive_error(self):
        from mass_accuracy_calculator import ppm_error

        assert ppm_error(500.0, 500.001) > 0

    def test_ppm_negative_error(self):
        from mass_accuracy_calculator import ppm_error

        assert ppm_error(500.0, 499.999) < 0

    def test_ppm_known_value(self):
        from mass_accuracy_calculator import ppm_error

        ppm = ppm_error(1000.0, 1000.001)
        assert abs(ppm - 1.0) < 0.001
