"""Tests for van_krevelen_data_generator."""

import csv
import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestComputeRatios:
    def test_glucose(self):
        from van_krevelen_data_generator import compute_ratios

        result = compute_ratios("C6H12O6")
        assert result["C"] == 6
        assert result["H"] == 12
        assert result["O"] == 6
        assert result["hc_ratio"] == pytest.approx(2.0, abs=0.01)
        assert result["oc_ratio"] == pytest.approx(1.0, abs=0.01)

    def test_palmitic_acid(self):
        from van_krevelen_data_generator import compute_ratios

        result = compute_ratios("C16H32O2")
        assert result["hc_ratio"] == pytest.approx(2.0, abs=0.01)
        assert result["oc_ratio"] == pytest.approx(0.125, abs=0.01)

    def test_no_carbon_raises(self):
        from van_krevelen_data_generator import compute_ratios

        with pytest.raises(ValueError, match="no carbon"):
            compute_ratios("H2O")


class TestClassifyCompound:
    def test_lipid_region(self):
        from van_krevelen_data_generator import classify_compound

        assert classify_compound(2.0, 0.125) == "lipids"

    def test_carbohydrate_region(self):
        from van_krevelen_data_generator import classify_compound

        assert classify_compound(2.0, 1.0) == "carbohydrates"

    def test_amino_acid_region(self):
        from van_krevelen_data_generator import classify_compound

        assert classify_compound(1.5, 0.5) == "amino_acids"

    def test_nucleotide_region(self):
        from van_krevelen_data_generator import classify_compound

        assert classify_compound(1.2, 0.7) == "nucleotides"

    def test_unclassified(self):
        from van_krevelen_data_generator import classify_compound

        assert classify_compound(0.5, 0.1) == "unclassified"


class TestProcessFormulas:
    def test_with_classification(self):
        from van_krevelen_data_generator import process_formulas

        results = process_formulas(["C6H12O6", "C16H32O2"], classify=True)
        assert len(results) == 2
        assert "class" in results[0]
        assert results[0]["class"] == "carbohydrates"
        assert results[1]["class"] == "lipids"

    def test_without_classification(self):
        from van_krevelen_data_generator import process_formulas

        results = process_formulas(["C6H12O6"], classify=False)
        assert "class" not in results[0]


class TestCLIRoundTrip:
    def test_roundtrip(self):
        from van_krevelen_data_generator import process_formulas

        formulas = ["C6H12O6", "C16H32O2", "C3H7NO2"]
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False, newline="") as inf:
            writer = csv.DictWriter(inf, fieldnames=["formula"], delimiter="\t")
            writer.writeheader()
            for f in formulas:
                writer.writerow({"formula": f})
            input_path = inf.name

        output_path = input_path.replace(".tsv", "_out.tsv")
        try:
            results = process_formulas(formulas, classify=True)
            fieldnames = ["formula", "C", "H", "O", "hc_ratio", "oc_ratio", "class"]
            with open(output_path, "w", newline="") as fh:
                writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
                writer.writeheader()
                writer.writerows(results)

            with open(output_path, newline="") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                rows = list(reader)
            assert len(rows) == 3
            assert "class" in rows[0]
        finally:
            os.unlink(input_path)
            if os.path.exists(output_path):
                os.unlink(output_path)
