"""
Protein Coverage Calculator
============================
Map identified peptides to proteins and calculate sequence coverage.

Features
--------
- Map peptides to protein sequences from FASTA
- Calculate per-protein sequence coverage percentage
- Report covered and uncovered regions

Usage
-----
    python protein_coverage_calculator.py --fasta proteins.fasta --peptides identified.tsv --output coverage.tsv
"""

import csv
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def load_fasta(fasta_path: str) -> dict:
    """Load proteins from a FASTA file.

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file.

    Returns
    -------
    dict
        Mapping of accession to sequence string.
    """
    entries = []
    oms.FASTAFile().load(fasta_path, entries)
    return {e.identifier: e.sequence for e in entries}


def calculate_coverage(proteins: dict, peptides: list) -> list:
    """Calculate sequence coverage for each protein.

    Parameters
    ----------
    proteins : dict
        Mapping of protein accession to sequence.
    peptides : list
        List of peptide sequence strings.

    Returns
    -------
    list
        List of dicts with coverage information per protein.
    """
    results = []
    for accession, protein_seq in proteins.items():
        prot_upper = protein_seq.upper()
        prot_len = len(prot_upper)
        covered = [False] * prot_len
        matched_peptides = []

        for pep in peptides:
            pep_upper = pep.strip().upper()
            start = 0
            while True:
                idx = prot_upper.find(pep_upper, start)
                if idx == -1:
                    break
                for i in range(idx, idx + len(pep_upper)):
                    covered[i] = True
                if pep not in matched_peptides:
                    matched_peptides.append(pep)
                start = idx + 1

        covered_count = sum(covered)
        coverage_pct = round(covered_count / prot_len * 100, 2) if prot_len > 0 else 0.0

        results.append({
            "accession": accession,
            "protein_length": prot_len,
            "covered_residues": covered_count,
            "coverage_percent": coverage_pct,
            "matched_peptides": len(matched_peptides),
            "peptides": ";".join(matched_peptides),
        })
    return results


@click.command(help="Calculate protein sequence coverage from peptides.")
@click.option("--fasta", required=True, help="Protein FASTA file.")
@click.option("--peptides", required=True, help="TSV file with 'sequence' column.")
@click.option("--output", default=None, help="Output TSV file.")
def main(fasta, peptides, output):
    """CLI entry point."""
    proteins = load_fasta(fasta)

    peptide_list = []
    with open(peptides) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            seq = row.get("sequence", "").strip()
            if seq:
                peptide_list.append(seq)

    results = calculate_coverage(proteins, peptide_list)

    if output:
        with open(output, "w", newline="") as fh:
            fieldnames = ["accession", "protein_length", "covered_residues", "coverage_percent",
                          "matched_peptides", "peptides"]
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(results)
        print(f"Results written to {output}")
    else:
        for r in results:
            print(f"{r['accession']}\t{r['coverage_percent']}%\t{r['matched_peptides']} peptides")


if __name__ == "__main__":
    main()
