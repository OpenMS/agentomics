"""
FragPipe Result Converter
=========================
Convert FragPipe psm.tsv to a standardized TSV format.

Maps FragPipe-specific column names to a common schema for downstream analysis.

Usage
-----
    python fragpipe_result_converter.py --input psm.tsv --output standardized.tsv
"""

import argparse
import csv
import sys

try:
    import pyopenms as oms  # noqa: F401
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

# Mapping from FragPipe column names to standard column names
COLUMN_MAP = {
    "Peptide": "peptide",
    "Modified Peptide": "modified_peptide",
    "Charge": "charge",
    "Calculated Peptide Mass": "mass",
    "Calibrated Observed Mass": "observed_mass",
    "Observed M/Z": "mz",
    "Retention": "rt",
    "Protein": "protein",
    "Protein Description": "protein_description",
    "Gene": "gene",
    "Hyperscore": "score",
    "Expectation": "expect",
    "PeptideProphet Probability": "probability",
    "Intensity": "intensity",
    "Spectrum": "spectrum",
    "Spectrum File": "raw_file",
    "Is Unique": "is_unique",
    "Mapped Proteins": "mapped_proteins",
}

STANDARD_FIELDS = [
    "peptide", "modified_peptide", "charge", "mass", "observed_mass", "mz", "rt",
    "protein", "protein_description", "gene",
    "score", "expect", "probability",
    "intensity", "spectrum", "raw_file",
    "is_unique", "mapped_proteins", "source",
]


def convert_fragpipe_psm(filepath: str) -> list:
    """Convert FragPipe psm.tsv to standardized format.

    Parameters
    ----------
    filepath:
        Path to FragPipe psm.tsv.

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
            for fp_col, std_col in COLUMN_MAP.items():
                std_row[std_col] = row.get(fp_col, "")
            std_row["source"] = "FragPipe"
            rows.append(std_row)
    return rows


def write_standardized(filepath: str, rows: list) -> None:
    """Write standardized results to TSV."""
    with open(filepath, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=STANDARD_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Convert FragPipe psm.tsv to standardized TSV.")
    parser.add_argument("--input", required=True, help="FragPipe psm.tsv file")
    parser.add_argument("--output", required=True, help="Output standardized TSV")
    args = parser.parse_args()

    rows = convert_fragpipe_psm(args.input)
    write_standardized(args.output, rows)

    print("Source: FragPipe")
    print(f"Total PSMs: {len(rows)}")
    print(f"Output written to {args.output}")


if __name__ == "__main__":
    main()
