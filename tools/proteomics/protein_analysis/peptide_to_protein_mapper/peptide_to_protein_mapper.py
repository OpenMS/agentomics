"""
Peptide to Protein Mapper
=========================
Map peptides to proteins by searching a FASTA database.

Uses pyopenms FASTAFile and AASequence to read FASTA entries and match
peptide sequences to protein sequences.

Usage
-----
    python peptide_to_protein_mapper.py --peptides peptides.tsv --fasta db.fasta --output mapped.tsv
"""

import csv

import click
import pyopenms as oms


def read_fasta(filepath: str) -> list:
    """Read a FASTA file using pyopenms.

    Parameters
    ----------
    filepath:
        Path to FASTA file.

    Returns
    -------
    list
        List of (accession, description, sequence) tuples.
    """
    entries = []
    fasta_entries = []
    fasta_file = oms.FASTAFile()
    fasta_file.load(filepath, fasta_entries)
    for entry in fasta_entries:
        entries.append((entry.identifier, entry.description, entry.sequence))
    return entries


def read_peptides(filepath: str) -> list:
    """Read a peptide list TSV.

    Expected columns: peptide (required), plus optional columns.

    Returns
    -------
    list
        List of dicts.
    """
    with open(filepath) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def map_peptides_to_proteins(peptides: list, fasta_entries: list) -> list:
    """Map each peptide to all proteins containing its sequence.

    Parameters
    ----------
    peptides:
        List of dicts with at least a 'peptide' key.
    fasta_entries:
        List of (accession, description, sequence) tuples.

    Returns
    -------
    list
        List of dicts with keys: peptide, protein, protein_description, start, end, is_unique.
    """
    results = []
    peptide_protein_counts = {}

    # First pass: count how many proteins each peptide maps to
    for pep_row in peptides:
        pep_seq = pep_row["peptide"].upper().strip()
        # Strip modifications for matching (simple bracket removal)
        clean_seq = _strip_modifications(pep_seq)
        count = 0
        for accession, description, prot_seq in fasta_entries:
            if clean_seq in prot_seq:
                count += 1
        peptide_protein_counts[pep_seq] = count

    # Second pass: build mappings
    for pep_row in peptides:
        pep_seq = pep_row["peptide"].upper().strip()
        clean_seq = _strip_modifications(pep_seq)
        is_unique = peptide_protein_counts.get(pep_seq, 0) == 1
        matched = False

        for accession, description, prot_seq in fasta_entries:
            start = prot_seq.find(clean_seq)
            if start >= 0:
                results.append({
                    "peptide": pep_seq,
                    "protein": accession,
                    "protein_description": description,
                    "start": start + 1,  # 1-based
                    "end": start + len(clean_seq),
                    "is_unique": is_unique,
                })
                matched = True

        if not matched:
            results.append({
                "peptide": pep_seq,
                "protein": "",
                "protein_description": "",
                "start": 0,
                "end": 0,
                "is_unique": False,
            })

    return results


def _strip_modifications(sequence: str) -> str:
    """Remove bracket-enclosed modifications from a peptide sequence."""
    result = []
    in_bracket = False
    for ch in sequence:
        if ch == "[":
            in_bracket = True
        elif ch == "]":
            in_bracket = False
        elif not in_bracket:
            result.append(ch)
    return "".join(result)


@click.command(help="Map peptides to proteins in a FASTA database.")
@click.option("--peptides", required=True, help="Input peptide TSV (must have 'peptide' column)")
@click.option("--fasta", required=True, help="FASTA database file")
@click.option("--output", required=True, help="Output TSV file")
def main(peptides, fasta, output):
    peptide_rows = read_peptides(peptides)
    fasta_entries = read_fasta(fasta)
    mappings = map_peptides_to_proteins(peptide_rows, fasta_entries)

    fieldnames = ["peptide", "protein", "protein_description", "start", "end", "is_unique"]
    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(mappings)

    n_mapped = sum(1 for m in mappings if m["protein"])
    n_unique = sum(1 for m in mappings if m["is_unique"])
    print(f"Peptides: {len(peptide_rows)}")
    print(f"Proteins in FASTA: {len(fasta_entries)}")
    print(f"Mappings: {n_mapped}")
    print(f"Unique peptides: {n_unique}")
    print(f"Output written to {output}")


if __name__ == "__main__":
    main()
