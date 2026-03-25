"""
Inclusion List Generator
========================
Generate instrument inclusion lists from peptide data for targeted MS experiments.

Supports output formats for Thermo and generic CSV. Calculates m/z values for
specified charge states using pyopenms AASequence.

Usage
-----
    python inclusion_list_generator.py --input peptides.tsv --format thermo --charge 2,3 --output inclusion.csv
"""

import argparse
import csv
import sys

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

PROTON = 1.007276


def read_peptides(filepath: str) -> list:
    """Read a peptide list TSV.

    Expected columns: peptide (required), plus optional: rt_start, rt_end, protein.

    Returns
    -------
    list
        List of dicts.
    """
    with open(filepath) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def calculate_mz(sequence: str, charge: int) -> float:
    """Calculate m/z for a peptide at the given charge state.

    Parameters
    ----------
    sequence:
        Peptide sequence (may include bracket modifications).
    charge:
        Charge state.

    Returns
    -------
    float
        Monoisotopic m/z value.
    """
    aa_seq = oms.AASequence.fromString(sequence)
    mono = aa_seq.getMonoWeight()
    return (mono + charge * PROTON) / charge


def generate_inclusion_list(
    peptides: list, charges: list, output_format: str = "thermo"
) -> list:
    """Generate an inclusion list for targeted MS.

    Parameters
    ----------
    peptides:
        List of dicts with at least a 'peptide' key.
    charges:
        List of charge states to include.
    output_format:
        'thermo' or 'generic'.

    Returns
    -------
    list
        List of dicts representing inclusion list entries.
    """
    output_format = output_format.lower()
    if output_format not in ("thermo", "generic"):
        raise ValueError(f"Unknown format: '{output_format}'. Choose 'thermo' or 'generic'.")

    entries = []
    for pep_row in peptides:
        sequence = pep_row["peptide"].strip()
        rt_start = pep_row.get("rt_start", "")
        rt_end = pep_row.get("rt_end", "")
        protein = pep_row.get("protein", "")

        for charge in charges:
            mz = calculate_mz(sequence, charge)

            if output_format == "thermo":
                entry = {
                    "Compound": sequence,
                    "Formula": "",
                    "Adduct": "",
                    "m/z": f"{mz:.4f}",
                    "z": str(charge),
                    "t start (min)": rt_start,
                    "t stop (min)": rt_end,
                    "Comment": protein,
                }
            else:
                entry = {
                    "peptide": sequence,
                    "mz": f"{mz:.4f}",
                    "charge": str(charge),
                    "rt_start": rt_start,
                    "rt_end": rt_end,
                    "protein": protein,
                }
            entries.append(entry)

    return entries


def main():
    parser = argparse.ArgumentParser(description="Generate instrument inclusion lists from peptide data.")
    parser.add_argument("--input", required=True, help="Input peptide TSV")
    parser.add_argument("--format", default="thermo", choices=["thermo", "generic"],
                        help="Output format (default: thermo)")
    parser.add_argument("--charge", default="2,3", help="Comma-separated charge states (default: 2,3)")
    parser.add_argument("--output", required=True, help="Output CSV file")
    args = parser.parse_args()

    charges = [int(c.strip()) for c in args.charge.split(",")]
    peptides = read_peptides(args.input)
    entries = generate_inclusion_list(peptides, charges, output_format=args.format)

    if not entries:
        print("No entries generated.")
        return

    fieldnames = list(entries[0].keys())
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(entries)

    print(f"Format: {args.format}")
    print(f"Charge states: {charges}")
    print(f"Peptides: {len(peptides)}")
    print(f"Entries: {len(entries)}")
    print(f"Output written to {args.output}")


if __name__ == "__main__":
    main()
