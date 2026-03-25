"""Tests for drug_metabolite_screener."""

import os
import tempfile

from conftest import requires_pyopenms


@requires_pyopenms
class TestPredictMetabolites:
    def test_phase1_oxidation(self):
        from drug_metabolite_screener import predict_metabolites

        results = predict_metabolites("C6H12O6", ["phase1"])
        reactions = [r["reaction"] for r in results]
        assert "oxidation" in reactions
        ox = next(r for r in results if r["reaction"] == "oxidation")
        assert ox["mass_shift"] > 0  # Adding O increases mass

    def test_phase2_glucuronidation(self):
        from drug_metabolite_screener import predict_metabolites

        results = predict_metabolites("C6H12O6", ["phase2"])
        reactions = [r["reaction"] for r in results]
        assert "glucuronidation" in reactions

    def test_combined_phases(self):
        from drug_metabolite_screener import predict_metabolites

        results = predict_metabolites("C17H14ClN3O", ["phase1", "phase2"])
        assert len(results) > 5  # Should have multiple reactions


@requires_pyopenms
class TestGetReactionTable:
    def test_phase1_only(self):
        from drug_metabolite_screener import get_reaction_table

        table = get_reaction_table(["phase1"])
        assert "oxidation" in table
        assert "glucuronidation" not in table

    def test_phase2_only(self):
        from drug_metabolite_screener import get_reaction_table

        table = get_reaction_table(["phase2"])
        assert "glucuronidation" in table
        assert "oxidation" not in table


@requires_pyopenms
class TestScreenMzML:
    def test_screen_synthetic_mzml(self):
        import pyopenms as oms
        from drug_metabolite_screener import predict_metabolites, screen_mzml

        metabolites = predict_metabolites("C6H12O6", ["phase1"])
        target_mass = metabolites[0]["exact_mass"]

        # Create synthetic mzML with a matching peak
        exp = oms.MSExperiment()
        spec = oms.MSSpectrum()
        spec.setMSLevel(1)
        spec.setRT(60.0)
        spec.set_peaks(([target_mass], [1000.0]))
        exp.addSpectrum(spec)

        with tempfile.NamedTemporaryFile(suffix=".mzML", delete=False) as tmp:
            oms.MzMLFile().store(tmp.name, exp)
            tmp_path = tmp.name

        try:
            matches = screen_mzml(tmp_path, metabolites, ppm=5.0)
            assert len(matches) >= 1
            assert abs(matches[0]["ppm_error"]) < 5.0
        finally:
            os.unlink(tmp_path)
