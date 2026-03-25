"""Tests for lipid_ecn_rt_predictor."""

import csv
import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestLipidEcnRtPredictor:
    def test_compute_ecn(self):
        from lipid_ecn_rt_predictor import compute_ecn

        assert compute_ecn(36, 2) == 32
        assert compute_ecn(34, 1) == 32
        assert compute_ecn(18, 0) == 18
        assert compute_ecn(20, 4) == 12

    def test_build_calibration_models(self):
        from lipid_ecn_rt_predictor import build_calibration_models

        standards = [
            {"lipid_class": "PC", "total_carbons": 32, "double_bonds": 0, "rt": 10.0},
            {"lipid_class": "PC", "total_carbons": 34, "double_bonds": 0, "rt": 12.0},
            {"lipid_class": "PC", "total_carbons": 36, "double_bonds": 0, "rt": 14.0},
        ]
        models = build_calibration_models(standards)
        assert "PC" in models
        assert abs(models["PC"]["r_value"] - 1.0) < 1e-6  # Perfect linear fit

    def test_build_calibration_single_point_skipped(self):
        from lipid_ecn_rt_predictor import build_calibration_models

        standards = [
            {"lipid_class": "PE", "total_carbons": 32, "double_bonds": 0, "rt": 10.0},
        ]
        models = build_calibration_models(standards)
        assert "PE" not in models

    def test_predict_rt(self):
        from lipid_ecn_rt_predictor import build_calibration_models, predict_rt

        standards = [
            {"lipid_class": "PC", "total_carbons": 32, "double_bonds": 0, "rt": 10.0},
            {"lipid_class": "PC", "total_carbons": 34, "double_bonds": 0, "rt": 12.0},
            {"lipid_class": "PC", "total_carbons": 36, "double_bonds": 0, "rt": 14.0},
        ]
        models = build_calibration_models(standards)

        lipids = [
            {"lipid_class": "PC", "total_carbons": 34, "double_bonds": 1},
        ]
        results = predict_rt(lipids, models)
        assert len(results) == 1
        assert results[0]["predicted_rt"] != "N/A"
        # ECN=32 -> RT=10.0 from the model (slope=1.0 per ECN)
        assert isinstance(results[0]["predicted_rt"], float)

    def test_predict_rt_missing_class(self):
        from lipid_ecn_rt_predictor import predict_rt

        lipids = [{"lipid_class": "SM", "total_carbons": 34, "double_bonds": 1}]
        results = predict_rt(lipids, {})
        assert results[0]["predicted_rt"] == "N/A"

    def test_full_pipeline(self):
        from lipid_ecn_rt_predictor import build_calibration_models, load_tsv, predict_rt, write_predictions

        with tempfile.TemporaryDirectory() as tmpdir:
            cal_path = os.path.join(tmpdir, "standards.tsv")
            with open(cal_path, "w", newline="") as fh:
                w = csv.DictWriter(fh, ["lipid_class", "total_carbons", "double_bonds", "rt"], delimiter="\t")
                w.writeheader()
                w.writerow({"lipid_class": "PC", "total_carbons": "32", "double_bonds": "0", "rt": "10.0"})
                w.writerow({"lipid_class": "PC", "total_carbons": "34", "double_bonds": "0", "rt": "12.0"})
                w.writerow({"lipid_class": "PC", "total_carbons": "36", "double_bonds": "0", "rt": "14.0"})

            lip_path = os.path.join(tmpdir, "lipids.tsv")
            with open(lip_path, "w", newline="") as fh:
                w = csv.DictWriter(fh, ["lipid_class", "total_carbons", "double_bonds"], delimiter="\t")
                w.writeheader()
                w.writerow({"lipid_class": "PC", "total_carbons": "34", "double_bonds": "1"})

            standards = load_tsv(cal_path)
            lipids = load_tsv(lip_path)
            models = build_calibration_models(standards)
            predictions = predict_rt(lipids, models)

            out_path = os.path.join(tmpdir, "predictions.tsv")
            write_predictions(predictions, out_path)
            assert os.path.exists(out_path)
