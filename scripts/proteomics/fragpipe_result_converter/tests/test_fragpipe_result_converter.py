"""Tests for fragpipe_result_converter."""

import csv

from conftest import requires_pyopenms
from fragpipe_result_converter import STANDARD_FIELDS, convert_fragpipe_psm, write_standardized


@requires_pyopenms
class TestFragpipeResultConverter:
    def _write_psm(self, tmp_path, rows):
        filepath = str(tmp_path / "psm.tsv")
        fieldnames = [
            "Peptide", "Modified Peptide", "Charge",
            "Calculated Peptide Mass", "Calibrated Observed Mass",
            "Observed M/Z", "Retention", "Protein", "Protein Description",
            "Gene", "Hyperscore", "Expectation",
            "PeptideProphet Probability", "Intensity",
            "Spectrum", "Spectrum File", "Is Unique", "Mapped Proteins",
        ]
        with open(filepath, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        return filepath

    def test_basic_conversion(self, tmp_path):
        filepath = self._write_psm(tmp_path, [
            {
                "Peptide": "PEPTIDEK", "Modified Peptide": "PEPTIDEK",
                "Charge": "2", "Calculated Peptide Mass": "900.0",
                "Calibrated Observed Mass": "900.1",
                "Observed M/Z": "450.5", "Retention": "25.5",
                "Protein": "sp|P12345|PROT1", "Protein Description": "Test",
                "Gene": "GEN1", "Hyperscore": "35.5",
                "Expectation": "1e-5",
                "PeptideProphet Probability": "0.999",
                "Intensity": "1e6", "Spectrum": "scan.1234.1234.2",
                "Spectrum File": "run1.mzML", "Is Unique": "true",
                "Mapped Proteins": "sp|P12345|PROT1",
            }
        ])
        rows = convert_fragpipe_psm(filepath)
        assert len(rows) == 1
        assert rows[0]["peptide"] == "PEPTIDEK"
        assert rows[0]["charge"] == "2"
        assert rows[0]["source"] == "FragPipe"

    def test_multiple_rows(self, tmp_path):
        filepath = self._write_psm(tmp_path, [
            {"Peptide": "PEP1", "Modified Peptide": "PEP1", "Charge": "2",
             "Calculated Peptide Mass": "800", "Calibrated Observed Mass": "800",
             "Observed M/Z": "400", "Retention": "20",
             "Protein": "P1", "Protein Description": "Prot1",
             "Gene": "G1", "Hyperscore": "30", "Expectation": "1e-4",
             "PeptideProphet Probability": "0.99", "Intensity": "1e6",
             "Spectrum": "scan.1", "Spectrum File": "run1.mzML",
             "Is Unique": "true", "Mapped Proteins": "P1"},
            {"Peptide": "PEP2", "Modified Peptide": "PEP2", "Charge": "3",
             "Calculated Peptide Mass": "900", "Calibrated Observed Mass": "900",
             "Observed M/Z": "300", "Retention": "30",
             "Protein": "P2", "Protein Description": "Prot2",
             "Gene": "G2", "Hyperscore": "25", "Expectation": "1e-3",
             "PeptideProphet Probability": "0.95", "Intensity": "5e5",
             "Spectrum": "scan.2", "Spectrum File": "run1.mzML",
             "Is Unique": "false", "Mapped Proteins": "P2;P3"},
        ])
        rows = convert_fragpipe_psm(filepath)
        assert len(rows) == 2

    def test_write_standardized(self, tmp_path):
        rows = [{"peptide": "PEPTIDEK", "charge": "2", "source": "FragPipe"}]
        outfile = str(tmp_path / "out.tsv")
        write_standardized(outfile, rows)
        with open(outfile) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            result = list(reader)
        assert len(result) == 1
        assert result[0]["source"] == "FragPipe"

    def test_standard_fields(self):
        assert "peptide" in STANDARD_FIELDS
        assert "source" in STANDARD_FIELDS
        assert "score" in STANDARD_FIELDS

    def test_missing_columns_handled(self, tmp_path):
        filepath = str(tmp_path / "minimal.tsv")
        with open(filepath, "w") as fh:
            fh.write("Peptide\tCharge\n")
            fh.write("PEPTIDEK\t2\n")
        rows = convert_fragpipe_psm(filepath)
        assert rows[0]["peptide"] == "PEPTIDEK"
        assert rows[0]["rt"] == ""
