"""Tests for hdx_back_exchange_estimator."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestHdxBackExchangeEstimator:
    def test_count_exchangeable_amides(self):
        from hdx_back_exchange_estimator import count_exchangeable_amides
        # AAALAAAK: 8 residues, 0 prolines -> 8-0-2 = 6
        assert count_exchangeable_amides("AAALAAAK") == 6

    def test_count_exchangeable_amides_with_proline(self):
        from hdx_back_exchange_estimator import count_exchangeable_amides
        # PPAAAA: 6 residues, 2 prolines -> 6-2-2 = 2
        assert count_exchangeable_amides("PPAAAA") == 2

    def test_get_peptide_mass(self):
        from hdx_back_exchange_estimator import get_peptide_mass
        mass = get_peptide_mass("PEPTIDEK")
        assert mass > 0

    def test_compute_theoretical_max(self):
        from hdx_back_exchange_estimator import (
            DEUTERIUM_MASS_SHIFT,
            compute_theoretical_max_deuterated_mass,
            count_exchangeable_amides,
            get_peptide_mass,
        )
        seq = "AAALAAAK"
        theoretical = compute_theoretical_max_deuterated_mass(seq)
        expected = get_peptide_mass(seq) + count_exchangeable_amides(seq) * DEUTERIUM_MASS_SHIFT
        assert abs(theoretical - expected) < 1e-6

    def test_compute_back_exchange_zero(self):
        from hdx_back_exchange_estimator import DEUTERIUM_MASS_SHIFT, compute_back_exchange, count_exchangeable_amides
        seq = "AAALAAAK"
        exchangeable = count_exchangeable_amides(seq)
        undeut = 700.0
        # Fully deuterated shows full shift -> 0% back-exchange
        fd = undeut + exchangeable * DEUTERIUM_MASS_SHIFT
        result = compute_back_exchange(seq, undeut, fd)
        assert abs(result["back_exchange_pct"]) < 0.01

    def test_compute_back_exchange_partial(self):
        from hdx_back_exchange_estimator import DEUTERIUM_MASS_SHIFT, compute_back_exchange, count_exchangeable_amides
        seq = "AAALAAAK"
        exchangeable = count_exchangeable_amides(seq)
        undeut = 700.0
        # Only 80% of theoretical shift observed -> 20% back-exchange
        fd = undeut + exchangeable * DEUTERIUM_MASS_SHIFT * 0.8
        result = compute_back_exchange(seq, undeut, fd)
        assert abs(result["back_exchange_pct"] - 20.0) < 0.1

    def test_flag_high_back_exchange(self):
        from hdx_back_exchange_estimator import flag_high_back_exchange
        results = [
            {"sequence": "AAA", "back_exchange_pct": 10.0},
            {"sequence": "BBB", "back_exchange_pct": 50.0},
        ]
        flagged = flag_high_back_exchange(results, max_backexchange=40.0)
        assert flagged[0]["exceeds_threshold"] == "NO"
        assert flagged[1]["exceeds_threshold"] == "YES"

    def test_write_output(self):
        from hdx_back_exchange_estimator import write_output
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "report.tsv")
            results = [{"sequence": "PEPTIDEK", "back_exchange_pct": 15.0, "exceeds_threshold": "NO"}]
            write_output(output_path, results)
            assert os.path.exists(output_path)

    def test_write_output_empty(self):
        from hdx_back_exchange_estimator import write_output
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "report.tsv")
            write_output(output_path, [])
            assert not os.path.exists(output_path)
