"""
Immunopeptide Filter
====================
Filter peptides for MHC class I or II binding by length and optional motif.

MHC-I peptides are typically 8-11 amino acids.
MHC-II peptides are typically 13-25 amino acids.

Usage
-----
    python immunopeptide_filter.py --input peptides.tsv --class-i --length-range 8-11 --output immunopeptides.tsv
    python immunopeptide_filter.py --input peptides.tsv --class-ii --output immunopeptides.tsv
"""

import argparse
import csv
import re
import sys
from typing import List, Optional, Tuple

try:
    import pyopenms as oms
except ImportError:
    sys.exit(
        "pyopenms is required. Install it with:  pip install pyopenms"
    )


def parse_length_range(range_str: str) -> Tuple[int, int]:
    """Parse a length range string like '8-11'.

    Parameters
    ----------
    range_str:
        Range string in format 'min-max'.

    Returns
    -------
    tuple
        (min_length, max_length)
    """
    parts = range_str.split("-")
    if len(parts) != 2:
        raise ValueError(f"Invalid range format: '{range_str}'. Expected 'min-max'.")
    return int(parts[0]), int(parts[1])


def filter_peptides(
    peptides: List[str],
    min_length: int = 8,
    max_length: int = 11,
    motif_pattern: Optional[str] = None,
) -> List[dict]:
    """Filter peptides by length and optional regex motif.

    Parameters
    ----------
    peptides:
        List of peptide sequences.
    min_length:
        Minimum peptide length (inclusive).
    max_length:
        Maximum peptide length (inclusive).
    motif_pattern:
        Optional regex pattern for motif filtering.

    Returns
    -------
    list
        List of dicts with sequence, length, and mass info for passing peptides.
    """
    results = []
    compiled_motif = re.compile(motif_pattern) if motif_pattern else None

    for seq_str in peptides:
        seq_str = seq_str.strip()
        if not seq_str:
            continue

        aa_seq = oms.AASequence.fromString(seq_str)
        length = aa_seq.size()

        if length < min_length or length > max_length:
            continue

        if compiled_motif and not compiled_motif.search(seq_str):
            continue

        results.append({
            "sequence": seq_str,
            "length": length,
            "monoisotopic_mass": aa_seq.getMonoWeight(),
        })

    return results


def read_peptides_from_tsv(input_path: str, column: str = "sequence") -> List[str]:
    """Read peptide sequences from a TSV file.

    Parameters
    ----------
    input_path:
        Path to input TSV file.
    column:
        Column name containing sequences.

    Returns
    -------
    list
        List of peptide sequence strings.
    """
    peptides = []
    with open(input_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if column in row and row[column].strip():
                peptides.append(row[column].strip())
    return peptides


def write_tsv(results: List[dict], output_path: str) -> None:
    """Write filtered results to TSV.

    Parameters
    ----------
    results:
        List of result dicts.
    output_path:
        Output file path.
    """
    if not results:
        with open(output_path, "w") as fh:
            fh.write("sequence\tlength\tmonoisotopic_mass\n")
        return
    fieldnames = ["sequence", "length", "monoisotopic_mass"]
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in results:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(
        description="Filter peptides for MHC-I/II by length and motif."
    )
    parser.add_argument("--input", required=True, help="Input TSV with peptide sequences")
    parser.add_argument("--column", default="sequence", help="Column name for sequences (default: sequence)")
    mhc_group = parser.add_mutually_exclusive_group()
    mhc_group.add_argument("--class-i", action="store_true", help="MHC class I defaults (8-11 aa)")
    mhc_group.add_argument("--class-ii", action="store_true", help="MHC class II defaults (13-25 aa)")
    parser.add_argument("--length-range", default=None, help="Custom length range, e.g. '8-11'")
    parser.add_argument("--motif", default=None, help="Regex motif pattern to filter by")
    parser.add_argument("--output", required=True, help="Output TSV file path")
    args = parser.parse_args()

    # Determine length range
    if args.length_range:
        min_len, max_len = parse_length_range(args.length_range)
    elif args.class_ii:
        min_len, max_len = 13, 25
    else:
        # Default to class I
        min_len, max_len = 8, 11

    peptides = read_peptides_from_tsv(args.input, column=args.column)
    print(f"Read {len(peptides)} peptides from {args.input}")

    results = filter_peptides(peptides, min_length=min_len, max_length=max_len, motif_pattern=args.motif)
    print(f"Passed filter: {len(results)} peptides (length {min_len}-{max_len})")

    write_tsv(results, args.output)
    print(f"Results written to {args.output}")


if __name__ == "__main__":
    main()
