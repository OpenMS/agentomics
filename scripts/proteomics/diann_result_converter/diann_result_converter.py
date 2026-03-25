"""
DIA-NN Result Converter
=======================
Convert DIA-NN report.tsv to a standardized TSV format.

Maps DIA-NN specific column names to a common schema for downstream analysis.

Usage
-----
    python diann_result_converter.py --input report.tsv --output standardized.tsv
"""

import argparse
import csv
import sys

try:
    import pyopenms as oms  # noqa: F401
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

# Mapping from DIA-NN column names to standard column names
COLUMN_MAP = {
    "Stripped.Sequence": "peptide",
    "Modified.Sequence": "modified_peptide",
    "Precursor.Charge": "charge",
    "Precursor.Mz": "mz",
    "RT": "rt",
    "Protein.Group": "protein",
    "Protein.Names": "protein_description",
    "Genes": "gene",
    "Q.Value": "qvalue",
    "PG.Q.Value": "pg_qvalue",
    "Global.Q.Value": "global_qvalue",
    "Precursor.Quantity": "intensity",
    "Run": "raw_file",
    "File.Name": "file_name",
}

STANDARD_FIELDS = [
    "peptide", "modified_peptide", "charge", "mz", "rt",
    "protein", "protein_description", "gene",
    "qvalue", "pg_qvalue", "global_qvalue",
    "intensity", "raw_file", "file_name", "source",
]


def convert_diann_report(filepath: str) -> list:
    """Convert DIA-NN report.tsv to standardized format.

    Parameters
    ----------
    filepath:
        Path to DIA-NN report.tsv.

    Returns
    -------
    list
        List of dicts with standardized column names.
    """
    rows = []
    with open(filepath) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            std_row = {}
            for diann_col, std_col in COLUMN_MAP.items():
                std_row[std_col] = row.get(diann_col, "")
            std_row["source"] = "DIA-NN"
            rows.append(std_row)
    return rows


def write_standardized(filepath: str, rows: list) -> None:
    """Write standardized results to TSV."""
    with open(filepath, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=STANDARD_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Convert DIA-NN report.tsv to standardized TSV.")
    parser.add_argument("--input", required=True, help="DIA-NN report.tsv file")
    parser.add_argument("--output", required=True, help="Output standardized TSV")
    args = parser.parse_args()

    rows = convert_diann_report(args.input)
    write_standardized(args.output, rows)

    print("Source: DIA-NN")
    print(f"Total precursors: {len(rows)}")
    print(f"Output written to {args.output}")


if __name__ == "__main__":
    main()
