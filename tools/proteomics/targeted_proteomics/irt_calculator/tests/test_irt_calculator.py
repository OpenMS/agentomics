"""Tests for irt_calculator."""


from conftest import requires_pyopenms


@requires_pyopenms
class TestIrtCalculator:
    def test_linear_regression_perfect(self):
        from irt_calculator import linear_regression

        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        y = [2.0, 4.0, 6.0, 8.0, 10.0]
        slope, intercept, r2 = linear_regression(x, y)
        assert abs(slope - 2.0) < 0.001
        assert abs(intercept - 0.0) < 0.001
        assert abs(r2 - 1.0) < 0.001

    def test_linear_regression_with_intercept(self):
        from irt_calculator import linear_regression

        x = [1.0, 2.0, 3.0]
        y = [3.0, 5.0, 7.0]  # y = 2x + 1
        slope, intercept, r2 = linear_regression(x, y)
        assert abs(slope - 2.0) < 0.001
        assert abs(intercept - 1.0) < 0.001

    def test_convert_to_irt(self):
        from irt_calculator import convert_to_irt

        irt = convert_to_irt(10.0, slope=2.0, intercept=5.0)
        assert abs(irt - 25.0) < 0.01

    def test_fit_irt_model(self):
        from irt_calculator import fit_irt_model

        reference = [
            {"sequence": "PEPTIDEK", "observed_rt": 10.0, "irt": 0.0},
            {"sequence": "ANOTHERPEPTIDE", "observed_rt": 20.0, "irt": 50.0},
            {"sequence": "THIRDPEPTIDE", "observed_rt": 30.0, "irt": 100.0},
        ]
        model = fit_irt_model(reference)
        assert model["r_squared"] > 0.99
        assert model["n_reference_peptides"] == 3

    def test_process_identifications(self):
        from irt_calculator import process_identifications

        model = {"slope": 5.0, "intercept": -50.0}
        idents = [
            {"sequence": "PEPTIDEK", "rt": 10.0},
            {"sequence": "ANOTHER", "rt": 20.0},
        ]
        results = process_identifications(idents, model)
        assert len(results) == 2
        assert results[0]["irt"] == 0.0  # 5*10 - 50
        assert results[1]["irt"] == 50.0  # 5*20 - 50

    def test_end_to_end_with_files(self):
        from irt_calculator import fit_irt_model, process_identifications

        reference = [
            {"sequence": "PEP1", "observed_rt": 5.0, "irt": -20.0},
            {"sequence": "PEP2", "observed_rt": 15.0, "irt": 30.0},
            {"sequence": "PEP3", "observed_rt": 25.0, "irt": 80.0},
        ]
        model = fit_irt_model(reference)
        idents = [{"sequence": "TEST", "rt": 10.0}]
        results = process_identifications(idents, model)
        assert len(results) == 1
        # Should be between -20 and 30
        assert -25 < results[0]["irt"] < 35
