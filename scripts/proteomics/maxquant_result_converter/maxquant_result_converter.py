"""
MaxQuant Result Converter
=========================
Convert MaxQuant evidence.txt to a standardized TSV format.

Maps MaxQuant-specific column names to a common schema suitable for
downstream analysis.

Usage
-----
    python maxquant_result_converter.py --input evidence.txt --output standardized.tsv
"""

import argparse
import csv
import sys

try:
    import pyopenms as oms  # noqa: F401
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

# Mapping from MaxQuant column names to standard column names
COLUMN_MAP = {
    "Sequence": "peptide",
    "Modified sequence": "modified_peptide",
    "Charge": "charge",
    "m/z": "mz",
    "Mass": "mass",
    "Retention time": "rt",
    "Proteins": "protein",
    "Leading razor protein": "leading_protein",
    "Gene names": "gene",
    "Score": "score",
    "PEP": "pep",
    "Intensity": "intensity",
    "Raw file": "raw_file",
    "Experiment": "experiment",
    "MS/MS scan number": "scan_number",
    "Reverse": "is_decoy",
    "Potential contaminant": "is_contaminant",
}

STANDARD_FIELDS = [
    "peptide", "modified_peptide", "charge", "mz", "mass", "rt",
    "protein", "leading_protein", "gene", "score", "pep",
    "intensity", "raw_file", "experiment", "scan_number",
    "is_decoy", "is_contaminant", "source",
]


def convert_maxquant_evidence(filepath: str) -> list:
    """Convert MaxQuant evidence.txt to standardized format.

    Parameters
    ----------
    filepath:
        Path to MaxQuant evidence.txt.

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
            for mq_col, std_col in COLUMN_MAP.items():
                value = row.get(mq_col, "")
                # Normalize decoy/contaminant flags
                if std_col == "is_decoy":
                    value = "true" if value == "+" else "false"
                elif std_col == "is_contaminant":
                    value = "true" if value == "+" else "false"
                std_row[std_col] = value
            std_row["source"] = "MaxQuant"
            rows.append(std_row)
    return rows


def write_standardized(filepath: str, rows: list) -> None:
    """Write standardized results to TSV."""
    with open(filepath, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=STANDARD_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Convert MaxQuant evidence.txt to standardized TSV.")
    parser.add_argument("--input", required=True, help="MaxQuant evidence.txt file")
    parser.add_argument("--output", required=True, help="Output standardized TSV")
    args = parser.parse_args()

    rows = convert_maxquant_evidence(args.input)
    write_standardized(args.output, rows)

    n_decoy = sum(1 for r in rows if r.get("is_decoy") == "true")
    print("Source: MaxQuant")
    print(f"Total PSMs: {len(rows)}")
    print(f"Decoy PSMs: {n_decoy}")
    print(f"Output written to {args.output}")


if __name__ == "__main__":
    main()
