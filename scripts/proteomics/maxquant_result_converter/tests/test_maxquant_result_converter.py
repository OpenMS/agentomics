"""Tests for maxquant_result_converter."""

import csv

from conftest import requires_pyopenms
from maxquant_result_converter import STANDARD_FIELDS, convert_maxquant_evidence, write_standardized


@requires_pyopenms
class TestMaxQuantResultConverter:
    def _write_evidence(self, tmp_path, rows):
        filepath = str(tmp_path / "evidence.txt")
        fieldnames = [
            "Sequence", "Modified sequence", "Charge", "m/z", "Mass",
            "Retention time", "Proteins", "Leading razor protein", "Gene names",
            "Score", "PEP", "Intensity", "Raw file", "Experiment",
            "MS/MS scan number", "Reverse", "Potential contaminant",
        ]
        with open(filepath, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        return filepath

    def test_basic_conversion(self, tmp_path):
        filepath = self._write_evidence(tmp_path, [
            {
                "Sequence": "PEPTIDEK", "Modified sequence": "_PEPTIDEK_",
                "Charge": "2", "m/z": "450.5", "Mass": "900.0",
                "Retention time": "25.5", "Proteins": "P12345",
                "Leading razor protein": "P12345", "Gene names": "GEN1",
                "Score": "120.5", "PEP": "0.001", "Intensity": "1e6",
                "Raw file": "run1", "Experiment": "exp1",
                "MS/MS scan number": "1234", "Reverse": "", "Potential contaminant": "",
            }
        ])
        rows = convert_maxquant_evidence(filepath)
        assert len(rows) == 1
        assert rows[0]["peptide"] == "PEPTIDEK"
        assert rows[0]["charge"] == "2"
        assert rows[0]["source"] == "MaxQuant"

    def test_decoy_flag(self, tmp_path):
        filepath = self._write_evidence(tmp_path, [
            {
                "Sequence": "PEPTIDEK", "Modified sequence": "_PEPTIDEK_",
                "Charge": "2", "m/z": "450.5", "Mass": "900.0",
                "Retention time": "25.5", "Proteins": "REV__P12345",
                "Leading razor protein": "REV__P12345", "Gene names": "",
                "Score": "50.0", "PEP": "0.5", "Intensity": "1e4",
                "Raw file": "run1", "Experiment": "exp1",
                "MS/MS scan number": "5678", "Reverse": "+", "Potential contaminant": "",
            }
        ])
        rows = convert_maxquant_evidence(filepath)
        assert rows[0]["is_decoy"] == "true"
        assert rows[0]["is_contaminant"] == "false"

    def test_contaminant_flag(self, tmp_path):
        filepath = self._write_evidence(tmp_path, [
            {
                "Sequence": "CONTAM", "Modified sequence": "_CONTAM_",
                "Charge": "1", "m/z": "300.0", "Mass": "299.0",
                "Retention time": "10.0", "Proteins": "CON__P00001",
                "Leading razor protein": "CON__P00001", "Gene names": "",
                "Score": "30.0", "PEP": "0.01", "Intensity": "5e5",
                "Raw file": "run1", "Experiment": "exp1",
                "MS/MS scan number": "999", "Reverse": "", "Potential contaminant": "+",
            }
        ])
        rows = convert_maxquant_evidence(filepath)
        assert rows[0]["is_contaminant"] == "true"

    def test_write_standardized(self, tmp_path):
        rows = [{"peptide": "PEPTIDEK", "charge": "2", "source": "MaxQuant"}]
        outfile = str(tmp_path / "out.tsv")
        write_standardized(outfile, rows)
        with open(outfile) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            result = list(reader)
        assert len(result) == 1
        assert result[0]["peptide"] == "PEPTIDEK"

    def test_standard_fields(self):
        assert "peptide" in STANDARD_FIELDS
        assert "source" in STANDARD_FIELDS

    def test_multiple_rows(self, tmp_path):
        filepath = self._write_evidence(tmp_path, [
            {"Sequence": "PEP1", "Modified sequence": "", "Charge": "2",
             "m/z": "400", "Mass": "800", "Retention time": "20",
             "Proteins": "P1", "Leading razor protein": "P1", "Gene names": "G1",
             "Score": "100", "PEP": "0.01", "Intensity": "1e6",
             "Raw file": "run1", "Experiment": "exp1",
             "MS/MS scan number": "1", "Reverse": "", "Potential contaminant": ""},
            {"Sequence": "PEP2", "Modified sequence": "", "Charge": "3",
             "m/z": "300", "Mass": "900", "Retention time": "30",
             "Proteins": "P2", "Leading razor protein": "P2", "Gene names": "G2",
             "Score": "90", "PEP": "0.02", "Intensity": "5e5",
             "Raw file": "run1", "Experiment": "exp1",
             "MS/MS scan number": "2", "Reverse": "", "Potential contaminant": ""},
        ])
        rows = convert_maxquant_evidence(filepath)
        assert len(rows) == 2
