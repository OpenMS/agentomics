"""Tests for suspect_screener."""

import csv
import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestSuspectScreener:
    def test_exact_mass_from_formula(self):
        from suspect_screener import exact_mass_from_formula

        mass = exact_mass_from_formula("C6H12O6")
        assert abs(mass - 180.0634) < 0.01

    def test_ppm_error_zero(self):
        from suspect_screener import ppm_error

        err = ppm_error(180.0634, 180.0634)
        assert abs(err) < 1e-6

    def test_ppm_error_positive(self):
        from suspect_screener import ppm_error

        err = ppm_error(180.0644, 180.0634)
        assert err > 0

    def test_ppm_error_negative(self):
        from suspect_screener import ppm_error

        err = ppm_error(180.0624, 180.0634)
        assert err < 0

    def test_screen_suspects_match(self):
        from suspect_screener import screen_suspects

        features = [{"feature_id": "F1", "mz": 180.0634, "rt": 120.0, "intensity": 1000.0}]
        suspects = [{"name": "Glucose", "formula": "C6H12O6", "exact_mass": 180.0634}]
        matches = screen_suspects(features, suspects, ppm_tolerance=5.0)
        assert len(matches) == 1
        assert matches[0]["suspect_name"] == "Glucose"
        assert abs(matches[0]["ppm_error"]) < 1.0

    def test_screen_suspects_no_match(self):
        from suspect_screener import screen_suspects

        features = [{"feature_id": "F1", "mz": 200.0, "rt": 120.0, "intensity": 1000.0}]
        suspects = [{"name": "Glucose", "formula": "C6H12O6", "exact_mass": 180.0634}]
        matches = screen_suspects(features, suspects, ppm_tolerance=5.0)
        assert len(matches) == 0

    def test_screen_suspects_sorted_by_abs_error(self):
        from suspect_screener import screen_suspects

        features = [{"feature_id": "F1", "mz": 180.0640, "rt": 120.0, "intensity": 1000.0}]
        suspects = [
            {"name": "A", "formula": "", "exact_mass": 180.0640},
            {"name": "B", "formula": "", "exact_mass": 180.0635},
        ]
        matches = screen_suspects(features, suspects, ppm_tolerance=10.0)
        assert len(matches) == 2
        assert matches[0]["abs_ppm_error"] <= matches[1]["abs_ppm_error"]

    def test_load_and_write_roundtrip(self):
        from suspect_screener import load_features, load_suspects, screen_suspects, write_matches

        with tempfile.TemporaryDirectory() as tmpdir:
            feat_path = os.path.join(tmpdir, "features.tsv")
            with open(feat_path, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=["feature_id", "mz", "rt", "intensity"], delimiter="\t")
                w.writeheader()
                w.writerow({"feature_id": "F1", "mz": "180.0634", "rt": "120", "intensity": "5000"})

            susp_path = os.path.join(tmpdir, "suspects.csv")
            with open(susp_path, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=["name", "formula", "exact_mass"])
                w.writeheader()
                w.writerow({"name": "Glucose", "formula": "C6H12O6", "exact_mass": "180.0634"})

            features = load_features(feat_path)
            suspects = load_suspects(susp_path)
            matches = screen_suspects(features, suspects, ppm_tolerance=5.0)

            out_path = os.path.join(tmpdir, "matches.tsv")
            write_matches(matches, out_path)
            assert os.path.exists(out_path)
            with open(out_path) as fh:
                lines = fh.readlines()
            assert len(lines) >= 2  # header + 1 match

    def test_load_suspects_computes_mass_from_formula(self):
        from suspect_screener import load_suspects

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "suspects.csv")
            with open(path, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=["name", "formula", "exact_mass"])
                w.writeheader()
                w.writerow({"name": "Glucose", "formula": "C6H12O6", "exact_mass": ""})

            suspects = load_suspects(path)
            assert len(suspects) == 1
            assert abs(suspects[0]["exact_mass"] - 180.0634) < 0.01
