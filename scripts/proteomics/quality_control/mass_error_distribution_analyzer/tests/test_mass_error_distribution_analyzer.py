"""Tests for mass_error_distribution_analyzer."""

from conftest import requires_pyopenms


@requires_pyopenms
class TestMassErrorDistributionAnalyzer:
    def test_compute_mass_errors(self):
        import pyopenms as oms
        from mass_error_distribution_analyzer import PROTON, compute_mass_errors

        aa = oms.AASequence.fromString("PEPTIDEK")
        theo_mass = aa.getMonoWeight()
        charge = 2
        theo_mz = (theo_mass + charge * PROTON) / charge

        rows = [{"sequence": "PEPTIDEK", "charge": "2", "precursor_mz": str(theo_mz)}]
        exp = oms.MSExperiment()

        errors = compute_mass_errors(rows, exp)
        assert len(errors) == 1
        assert abs(errors[0]["error_ppm"]) < 0.01

    def test_mass_error_with_offset(self):
        import pyopenms as oms
        from mass_error_distribution_analyzer import PROTON, compute_mass_errors

        aa = oms.AASequence.fromString("PEPTIDEK")
        theo_mass = aa.getMonoWeight()
        charge = 2
        theo_mz = (theo_mass + charge * PROTON) / charge
        obs_mz = theo_mz + 0.001  # 1 mDa offset

        rows = [{"sequence": "PEPTIDEK", "charge": "2", "precursor_mz": str(obs_mz)}]
        exp = oms.MSExperiment()

        errors = compute_mass_errors(rows, exp)
        assert len(errors) == 1
        assert errors[0]["error_ppm"] > 0

    def test_summarize_errors(self):
        from mass_error_distribution_analyzer import summarize_errors

        errors = [
            {"error_ppm": 1.0, "error_da": 0.0005},
            {"error_ppm": -1.0, "error_da": -0.0005},
            {"error_ppm": 0.5, "error_da": 0.00025},
        ]
        summary = summarize_errors(errors)
        assert summary["count"] == 3
        assert abs(summary["ppm_mean"] - (1.0 - 1.0 + 0.5) / 3) < 0.01

    def test_empty_errors(self):
        from mass_error_distribution_analyzer import summarize_errors

        summary = summarize_errors([])
        assert summary["count"] == 0
