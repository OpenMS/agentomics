"""Tests for formula_validator_golden_rules."""


import pytest

pytest.importorskip("pyopenms")


class TestGetElementCounts:
    def test_glucose(self):
        from formula_validator_golden_rules import get_element_counts

        counts = get_element_counts("C6H12O6")
        assert counts["C"] == 6
        assert counts["H"] == 12
        assert counts["O"] == 6

    def test_with_nitrogen(self):
        from formula_validator_golden_rules import get_element_counts

        counts = get_element_counts("C3H7NO2")
        assert counts["N"] == 1


class TestComputeRDBE:
    def test_benzene(self):
        from formula_validator_golden_rules import compute_rdbe

        # C6H6: RDBE = (12 + 2 - 6) / 2 = 4
        assert compute_rdbe({"C": 6, "H": 6}) == 4.0

    def test_methane(self):
        from formula_validator_golden_rules import compute_rdbe

        # CH4: RDBE = (2 + 2 - 4) / 2 = 0
        assert compute_rdbe({"C": 1, "H": 4}) == 0.0

    def test_glucose(self):
        from formula_validator_golden_rules import compute_rdbe

        # C6H12O6: RDBE = (12 + 2 - 12) / 2 = 1
        assert compute_rdbe({"C": 6, "H": 12, "O": 6}) == 1.0


class TestIndividualRules:
    def test_rdbe_nonneg_pass(self):
        from formula_validator_golden_rules import check_rdbe_nonnegative

        assert check_rdbe_nonnegative({"C": 6, "H": 6}) is True

    def test_rdbe_nonneg_fail(self):
        from formula_validator_golden_rules import check_rdbe_nonnegative

        # Very high H count makes RDBE negative
        assert check_rdbe_nonnegative({"C": 1, "H": 10}) is False

    def test_hc_ratio_pass(self):
        from formula_validator_golden_rules import check_hc_ratio

        assert check_hc_ratio({"C": 6, "H": 12}) is True  # H/C = 2.0

    def test_hc_ratio_fail_too_low(self):
        from formula_validator_golden_rules import check_hc_ratio

        assert check_hc_ratio({"C": 10, "H": 1}) is False  # H/C = 0.1

    def test_nc_ratio_pass(self):
        from formula_validator_golden_rules import check_nc_ratio

        assert check_nc_ratio({"C": 3, "N": 1}) is True  # N/C = 0.33

    def test_oc_ratio_fail(self):
        from formula_validator_golden_rules import check_oc_ratio

        assert check_oc_ratio({"C": 1, "O": 5}) is False  # O/C = 5.0

    def test_sc_ratio_pass(self):
        from formula_validator_golden_rules import check_sc_ratio

        assert check_sc_ratio({"C": 10, "S": 1}) is True

    def test_pc_ratio_pass(self):
        from formula_validator_golden_rules import check_pc_ratio

        assert check_pc_ratio({"C": 10, "P": 1}) is True


class TestValidateFormula:
    def test_glucose_valid(self):
        from formula_validator_golden_rules import validate_formula

        result = validate_formula("C6H12O6", ["all"])
        assert result["valid"] is True

    def test_alanine_valid(self):
        from formula_validator_golden_rules import validate_formula

        result = validate_formula("C3H7NO2", ["all"])
        assert result["valid"] is True

    def test_specific_rules(self):
        from formula_validator_golden_rules import validate_formula

        result = validate_formula("C6H12O6", ["rdbe", "hc"])
        assert "rule_rdbe" in result
        assert "rule_hc" in result
        assert "rule_nc" not in result


class TestValidateFormulas:
    def test_batch(self):
        from formula_validator_golden_rules import validate_formulas

        results = validate_formulas(["C6H12O6", "C3H7NO2"], ["all"])
        assert len(results) == 2
        assert all(r["valid"] for r in results)
