"""
Peptide Uniqueness Checker
==========================
Check if peptides are proteotypic (unique to a single protein) in a FASTA database.

Features
--------
- Map peptides to all matching proteins in a FASTA file
- Flag proteotypic (unique) vs shared peptides
- Report all matching protein accessions

Usage
-----
    python peptide_uniqueness_checker.py --peptides peptides.tsv --fasta db.fasta --output uniqueness.tsv
"""

import argparse
import csv
import sys

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def load_fasta(fasta_path: str) -> list:
    """Load protein entries from a FASTA file.

    Parameters
    ----------
    fasta_path : str
        Path to FASTA file.

    Returns
    -------
    list
        List of (accession, sequence) tuples.
    """
    entries = []
    fasta_file = oms.FASTAFile()
    fasta_entries = []
    fasta_file.load(fasta_path, fasta_entries)
    for entry in fasta_entries:
        entries.append((entry.identifier, entry.sequence))
    return entries


def check_uniqueness(peptides: list, fasta_path: str) -> list:
    """Check peptide uniqueness against a FASTA database.

    Parameters
    ----------
    peptides : list
        List of peptide sequence strings.
    fasta_path : str
        Path to the FASTA file.

    Returns
    -------
    list
        List of dicts with keys: peptide, proteins, protein_count, is_proteotypic.
    """
    proteins = load_fasta(fasta_path)
    results = []
    for pep in peptides:
        pep_upper = pep.strip().upper()
        matching = []
        for accession, seq in proteins:
            if pep_upper in seq.upper():
                matching.append(accession)
        results.append({
            "peptide": pep,
            "proteins": ";".join(matching),
            "protein_count": len(matching),
            "is_proteotypic": len(matching) == 1,
        })
    return results


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(description="Check peptide uniqueness in a FASTA database.")
    parser.add_argument("--peptides", required=True, help="TSV file with 'sequence' column, or comma-separated list.")
    parser.add_argument("--fasta", required=True, help="FASTA database file.")
    parser.add_argument("--output", help="Output TSV file.")
    args = parser.parse_args()

    # Load peptides
    peptide_list = []
    try:
        with open(args.peptides) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                seq = row.get("sequence", "").strip()
                if seq:
                    peptide_list.append(seq)
    except (FileNotFoundError, KeyError):
        peptide_list = [p.strip() for p in args.peptides.split(",") if p.strip()]

    results = check_uniqueness(peptide_list, args.fasta)

    if args.output:
        with open(args.output, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=["peptide", "proteins", "protein_count", "is_proteotypic"],
                                    delimiter="\t")
            writer.writeheader()
            writer.writerows(results)
        print(f"Results written to {args.output}")
    else:
        for r in results:
            status = "proteotypic" if r["is_proteotypic"] else "shared"
            print(f"{r['peptide']}\t{status}\t{r['protein_count']}\t{r['proteins']}")


if __name__ == "__main__":
    main()
