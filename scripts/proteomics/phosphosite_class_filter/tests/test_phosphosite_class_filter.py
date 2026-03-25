"""Tests for phosphosite_class_filter."""

import csv
import os
import tempfile

from conftest import requires_pyopenms


@requires_pyopenms
class TestPhosphositeClassFilter:
    def test_classify_class1(self):
        from phosphosite_class_filter import classify_phosphosite
        assert classify_phosphosite(0.90, 0.75) == "Class I"

    def test_classify_class2(self):
        from phosphosite_class_filter import classify_phosphosite
        assert classify_phosphosite(0.60, 0.75) == "Class II"

    def test_classify_class3(self):
        from phosphosite_class_filter import classify_phosphosite
        assert classify_phosphosite(0.30, 0.75) == "Class III"

    def test_classify_boundary(self):
        from phosphosite_class_filter import classify_phosphosite
        assert classify_phosphosite(0.75, 0.75) == "Class I"
        assert classify_phosphosite(0.50, 0.75) == "Class II"

    def test_validate_phosphopeptide(self):
        from phosphosite_class_filter import validate_phosphopeptide
        assert validate_phosphopeptide("PEPTIDEK") is True

    def test_classify_phosphosites(self):
        from phosphosite_class_filter import classify_phosphosites
        rows = [
            {"peptide": "PEPTIDEK", "protein": "P1", "site": "S5", "localization_prob": "0.90",
             "modification": "Phospho"},
            {"peptide": "ACDEFGH", "protein": "P2", "site": "T3", "localization_prob": "0.60",
             "modification": "Phospho"},
            {"peptide": "KLMNPQR", "protein": "P3", "site": "Y7", "localization_prob": "0.30",
             "modification": "Phospho"},
        ]
        classified, summary = classify_phosphosites(rows, 0.75)
        assert len(classified) == 3
        assert summary["Class I"] == 1
        assert summary["Class II"] == 1
        assert summary["Class III"] == 1

    def test_enrichment_efficiency(self):
        from phosphosite_class_filter import compute_enrichment_efficiency
        summary = {"Class I": 3, "Class II": 1, "Class III": 1, "total": 5}
        assert abs(compute_enrichment_efficiency(summary) - 0.6) < 1e-6

    def test_enrichment_efficiency_empty(self):
        from phosphosite_class_filter import compute_enrichment_efficiency
        summary = {"Class I": 0, "Class II": 0, "Class III": 0, "total": 0}
        assert compute_enrichment_efficiency(summary) == 0.0

    def test_read_write_roundtrip(self):
        from phosphosite_class_filter import classify_phosphosites, read_input, write_output

        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "input.tsv")
            output_path = os.path.join(tmpdir, "output.tsv")
            with open(input_path, "w", newline="") as f:
                writer = csv.DictWriter(
                    f, fieldnames=["peptide", "protein", "site", "localization_prob", "modification"],
                    delimiter="\t"
                )
                writer.writeheader()
                writer.writerow({
                    "peptide": "PEPTIDEK", "protein": "P1", "site": "S5",
                    "localization_prob": "0.85", "modification": "Phospho"
                })

            rows = read_input(input_path)
            classified, summary = classify_phosphosites(rows)
            write_output(output_path, classified, summary)
            assert os.path.exists(output_path)

            with open(output_path) as f:
                reader = csv.DictReader(f, delimiter="\t")
                result_rows = list(reader)
            assert len(result_rows) == 1
            assert result_rows[0]["site_class"] == "Class I"
