"""Tests for charge_state_predictor."""

import pytest

pytest.importorskip("pyopenms")


class TestChargeStatePredictor:
    def test_basic_sites_counting(self):
        from charge_state_predictor import count_basic_sites

        sites = count_basic_sites("PEPTIDEK")
        assert sites["nterm"] == 1
        assert sites["K"] == 1
        assert sites["total"] == 2

    def test_arginine_sites(self):
        from charge_state_predictor import count_basic_sites

        sites = count_basic_sites("PEPTIDERR")
        assert sites["R"] == 2

    def test_histidine_sites(self):
        from charge_state_predictor import count_basic_sites

        sites = count_basic_sites("HPEPTIDEH")
        assert sites["H"] == 2

    def test_predict_charge_states(self):
        from charge_state_predictor import predict_charge_states

        result = predict_charge_states("PEPTIDEK", ph=2.0)
        assert len(result["charge_states"]) > 0
        assert result["monoisotopic_mass"] > 0
        # Probabilities should sum to ~1
        total_prob = sum(cs["probability"] for cs in result["charge_states"])
        assert abs(total_prob - 1.0) < 0.01

    def test_more_basic_residues_higher_charge(self):
        from charge_state_predictor import predict_charge_states

        few_basic = predict_charge_states("PEPTIDEA", ph=2.0)
        many_basic = predict_charge_states("KPEPTIDEKRKH", ph=2.0)
        assert many_basic["expected_charge"] > few_basic["expected_charge"]

    def test_low_ph_higher_charge(self):
        from charge_state_predictor import predict_charge_states

        low_ph = predict_charge_states("PEPTIDEK", ph=1.0)
        high_ph = predict_charge_states("PEPTIDEK", ph=7.0)
        assert low_ph["expected_charge"] >= high_ph["expected_charge"]

    def test_mz_values(self):
        from charge_state_predictor import PROTON, predict_charge_states

        result = predict_charge_states("PEPTIDEK", ph=2.0)
        for cs in result["charge_states"]:
            expected_mz = (result["monoisotopic_mass"] + cs["charge"] * PROTON) / cs["charge"]
            assert abs(cs["mz"] - expected_mz) < 0.001

    def test_most_likely_charge(self):
        from charge_state_predictor import predict_charge_states

        result = predict_charge_states("PEPTIDEK", ph=2.0)
        assert result["most_likely_charge"] >= 1
