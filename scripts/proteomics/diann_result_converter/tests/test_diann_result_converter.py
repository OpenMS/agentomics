"""Tests for diann_result_converter."""

import csv

from conftest import requires_pyopenms
from diann_result_converter import STANDARD_FIELDS, convert_diann_report, write_standardized


@requires_pyopenms
class TestDiannResultConverter:
    def _write_report(self, tmp_path, rows):
        filepath = str(tmp_path / "report.tsv")
        fieldnames = [
            "Stripped.Sequence", "Modified.Sequence", "Precursor.Charge",
            "Precursor.Mz", "RT", "Protein.Group", "Protein.Names",
            "Genes", "Q.Value", "PG.Q.Value", "Global.Q.Value",
            "Precursor.Quantity", "Run", "File.Name",
        ]
        with open(filepath, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        return filepath

    def test_basic_conversion(self, tmp_path):
        filepath = self._write_report(tmp_path, [
            {
                "Stripped.Sequence": "PEPTIDEK", "Modified.Sequence": "PEPTIDEK",
                "Precursor.Charge": "2", "Precursor.Mz": "450.5",
                "RT": "25.5", "Protein.Group": "P12345",
                "Protein.Names": "Protein 1", "Genes": "GEN1",
                "Q.Value": "0.001", "PG.Q.Value": "0.005",
                "Global.Q.Value": "0.01",
                "Precursor.Quantity": "1e6", "Run": "run1",
                "File.Name": "run1.mzML",
            }
        ])
        rows = convert_diann_report(filepath)
        assert len(rows) == 1
        assert rows[0]["peptide"] == "PEPTIDEK"
        assert rows[0]["charge"] == "2"
        assert rows[0]["source"] == "DIA-NN"

    def test_multiple_rows(self, tmp_path):
        filepath = self._write_report(tmp_path, [
            {"Stripped.Sequence": "PEP1", "Modified.Sequence": "PEP1",
             "Precursor.Charge": "2", "Precursor.Mz": "400",
             "RT": "20", "Protein.Group": "P1", "Protein.Names": "Prot1",
             "Genes": "G1", "Q.Value": "0.01", "PG.Q.Value": "0.02",
             "Global.Q.Value": "0.03", "Precursor.Quantity": "1e6",
             "Run": "run1", "File.Name": "run1.mzML"},
            {"Stripped.Sequence": "PEP2", "Modified.Sequence": "PEP2",
             "Precursor.Charge": "3", "Precursor.Mz": "300",
             "RT": "30", "Protein.Group": "P2", "Protein.Names": "Prot2",
             "Genes": "G2", "Q.Value": "0.02", "PG.Q.Value": "0.03",
             "Global.Q.Value": "0.04", "Precursor.Quantity": "5e5",
             "Run": "run1", "File.Name": "run1.mzML"},
        ])
        rows = convert_diann_report(filepath)
        assert len(rows) == 2

    def test_write_standardized(self, tmp_path):
        rows = [{"peptide": "PEPTIDEK", "charge": "2", "source": "DIA-NN"}]
        outfile = str(tmp_path / "out.tsv")
        write_standardized(outfile, rows)
        with open(outfile) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            result = list(reader)
        assert len(result) == 1
        assert result[0]["source"] == "DIA-NN"

    def test_standard_fields(self):
        assert "peptide" in STANDARD_FIELDS
        assert "source" in STANDARD_FIELDS
        assert "qvalue" in STANDARD_FIELDS

    def test_missing_columns_handled(self, tmp_path):
        filepath = str(tmp_path / "minimal.tsv")
        with open(filepath, "w") as fh:
            fh.write("Stripped.Sequence\tPrecursor.Charge\n")
            fh.write("PEPTIDEK\t2\n")
        rows = convert_diann_report(filepath)
        assert rows[0]["peptide"] == "PEPTIDEK"
        assert rows[0]["rt"] == ""  # Missing column => empty string
