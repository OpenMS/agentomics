"""Tests for metabolite_class_predictor."""

import csv
import os
import tempfile

import pytest
from conftest import requires_pyopenms


@requires_pyopenms
class TestMetaboliteClassPredictor:
    def test_get_element_counts(self):
        from metabolite_class_predictor import get_element_counts

        counts = get_element_counts("C6H12O6")
        assert counts["C"] == 6
        assert counts["H"] == 12
        assert counts["O"] == 6

    def test_get_element_counts_with_nitrogen(self):
        from metabolite_class_predictor import get_element_counts

        counts = get_element_counts("C2H5NO2")  # Glycine
        assert counts["C"] == 2
        assert counts["N"] == 1

    def test_compute_exact_mass(self):
        from metabolite_class_predictor import compute_exact_mass

        mass = compute_exact_mass("C6H12O6")
        assert abs(mass - 180.0634) < 0.01

    def test_compute_mass_defect(self):
        from metabolite_class_predictor import compute_mass_defect

        defect = compute_mass_defect(180.0634)
        assert 0.0 < defect < 1.0
        assert abs(defect - 0.0634) < 0.01

    def test_compute_rdbe(self):
        from metabolite_class_predictor import compute_rdbe

        # Benzene C6H6: RDBE = 1 + 6 - 6/2 = 4
        rdbe = compute_rdbe({"C": 6, "H": 6})
        assert abs(rdbe - 4.0) < 0.01

        # Glucose C6H12O6: RDBE = 1 + 6 - 12/2 = 1
        rdbe = compute_rdbe({"C": 6, "H": 12, "O": 6})
        assert abs(rdbe - 1.0) < 0.01

    def test_compute_element_ratios(self):
        from metabolite_class_predictor import compute_element_ratios

        ratios = compute_element_ratios({"C": 6, "H": 12, "O": 6})
        assert ratios["hc_ratio"] == pytest.approx(2.0)
        assert ratios["oc_ratio"] == pytest.approx(1.0)

    def test_compute_element_ratios_no_carbon(self):
        from metabolite_class_predictor import compute_element_ratios

        ratios = compute_element_ratios({"H": 2, "O": 1})
        assert ratios["hc_ratio"] is None
        assert ratios["oc_ratio"] is None

    def test_classify_carbohydrate(self):
        from metabolite_class_predictor import classify_metabolite

        # Glucose C6H12O6: H:C=2.0, O:C=1.0 -> Carbohydrate
        result = classify_metabolite("C6H12O6")
        assert result["predicted_class"] == "Carbohydrate"

    def test_classify_lipid(self):
        from metabolite_class_predictor import classify_metabolite

        # Palmitic acid C16H32O2: H:C=2.0, O:C=0.125 -> Lipid
        result = classify_metabolite("C16H32O2")
        assert result["predicted_class"] == "Lipid"

    def test_classify_amino_acid(self):
        from metabolite_class_predictor import classify_metabolite

        # Alanine C3H7NO2: H:C=2.33, O:C=0.67, has N
        result = classify_metabolite("C3H7NO2")
        assert "Amino acid" in result["predicted_class"] or "Peptide" in result["predicted_class"]

    def test_classify_returns_all_fields(self):
        from metabolite_class_predictor import classify_metabolite

        result = classify_metabolite("C6H12O6")
        expected_keys = [
            "formula", "exact_mass", "mass_defect", "rdbe",
            "hc_ratio", "oc_ratio", "has_nitrogen", "has_sulfur",
            "has_phosphorus", "predicted_class", "confidence",
        ]
        for key in expected_keys:
            assert key in result

    def test_classify_batch(self):
        from metabolite_class_predictor import classify_batch

        formulas = [{"formula": "C6H12O6"}, {"formula": "C16H32O2"}]
        results = classify_batch(formulas)
        assert len(results) == 2

    def test_classify_batch_skips_empty(self):
        from metabolite_class_predictor import classify_batch

        formulas = [{"formula": ""}, {"formula": "C6H12O6"}]
        results = classify_batch(formulas)
        assert len(results) == 1

    def test_full_pipeline(self):
        from metabolite_class_predictor import classify_batch, load_formulas, write_predictions

        with tempfile.TemporaryDirectory() as tmpdir:
            in_path = os.path.join(tmpdir, "formulas.tsv")
            with open(in_path, "w", newline="") as fh:
                w = csv.DictWriter(fh, ["formula"], delimiter="\t")
                w.writeheader()
                w.writerow({"formula": "C6H12O6"})
                w.writerow({"formula": "C16H32O2"})

            formulas = load_formulas(in_path)
            predictions = classify_batch(formulas)

            out_path = os.path.join(tmpdir, "predictions.tsv")
            write_predictions(predictions, out_path)
            assert os.path.exists(out_path)
