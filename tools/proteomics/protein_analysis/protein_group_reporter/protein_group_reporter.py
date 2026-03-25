"""
Protein Group Reporter
======================
Parse protein groups from peptide-level data and a FASTA database, then
report a clean protein group table with peptide counts and sequence coverage.

Usage
-----
    python protein_group_reporter.py --input peptides.tsv --fasta db.fasta --output groups.tsv
"""

import csv
from typing import Dict, List, Set

import click
import pyopenms as oms


def load_fasta(fasta_path: str) -> Dict[str, str]:
    """Load a FASTA file and return a dict mapping accession to sequence.

    Parameters
    ----------
    fasta_path:
        Path to a FASTA file.

    Returns
    -------
    dict
        Mapping of accession to protein sequence string.
    """
    entries = []
    fasta_file = oms.FASTAFile()
    fasta_file.load(fasta_path, entries)
    proteins = {}
    for entry in entries:
        proteins[entry.identifier] = entry.sequence
    return proteins


def map_peptides_to_proteins(
    peptides: List[str],
    proteins: Dict[str, str],
) -> Dict[str, Set[str]]:
    """Map peptides to proteins by substring matching.

    Parameters
    ----------
    peptides:
        List of peptide sequences.
    proteins:
        Dict mapping accession to protein sequence.

    Returns
    -------
    dict
        Mapping of protein accession to set of matched peptide sequences.
    """
    protein_peptides: Dict[str, Set[str]] = {}
    for pep in peptides:
        pep = pep.strip()
        if not pep:
            continue
        for acc, prot_seq in proteins.items():
            if pep in prot_seq:
                if acc not in protein_peptides:
                    protein_peptides[acc] = set()
                protein_peptides[acc].add(pep)
    return protein_peptides


def compute_sequence_coverage(protein_seq: str, peptides: Set[str]) -> float:
    """Compute sequence coverage of a protein by its peptides.

    Parameters
    ----------
    protein_seq:
        Full protein sequence.
    peptides:
        Set of peptide sequences mapped to this protein.

    Returns
    -------
    float
        Fraction of protein sequence covered (0.0 to 1.0).
    """
    if not protein_seq:
        return 0.0
    covered = [False] * len(protein_seq)
    for pep in peptides:
        start = 0
        while True:
            idx = protein_seq.find(pep, start)
            if idx == -1:
                break
            for i in range(idx, idx + len(pep)):
                covered[i] = True
            start = idx + 1
    return sum(covered) / len(protein_seq)


def build_protein_groups(
    protein_peptides: Dict[str, Set[str]],
    proteins: Dict[str, str],
) -> List[dict]:
    """Build protein group report.

    Groups proteins that share the exact same set of peptides.

    Parameters
    ----------
    protein_peptides:
        Mapping of accession to peptide set.
    proteins:
        Mapping of accession to protein sequence.

    Returns
    -------
    list
        List of protein group dicts.
    """
    # Group proteins with identical peptide sets
    peptide_set_to_accessions: Dict[frozenset, List[str]] = {}
    for acc, peps in protein_peptides.items():
        key = frozenset(peps)
        if key not in peptide_set_to_accessions:
            peptide_set_to_accessions[key] = []
        peptide_set_to_accessions[key].append(acc)

    results = []
    for pep_set, accessions in peptide_set_to_accessions.items():
        accessions_sorted = sorted(accessions)
        lead = accessions_sorted[0]
        coverage = compute_sequence_coverage(proteins.get(lead, ""), pep_set)
        results.append({
            "protein_group": ";".join(accessions_sorted),
            "lead_protein": lead,
            "n_proteins": len(accessions_sorted),
            "n_peptides": len(pep_set),
            "peptides": ";".join(sorted(pep_set)),
            "sequence_coverage": round(coverage, 4),
        })

    results.sort(key=lambda x: x["n_peptides"], reverse=True)
    return results


def read_peptides_from_tsv(input_path: str, column: str = "sequence") -> List[str]:
    """Read peptide sequences from a TSV file.

    Parameters
    ----------
    input_path:
        Path to input TSV.
    column:
        Column name for peptide sequences.

    Returns
    -------
    list
        List of peptide sequences.
    """
    peptides = []
    with open(input_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if column in row and row[column].strip():
                peptides.append(row[column].strip())
    return peptides


def write_tsv(results: List[dict], output_path: str) -> None:
    """Write protein group results to TSV.

    Parameters
    ----------
    results:
        List of protein group dicts.
    output_path:
        Output file path.
    """
    fieldnames = ["protein_group", "lead_protein", "n_proteins", "n_peptides", "peptides", "sequence_coverage"]
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in results:
            writer.writerow(row)


@click.command(help="Parse protein groups from peptide data and report clean table.")
@click.option("--input", "input", required=True, help="Input TSV with peptide sequences")
@click.option("--fasta", required=True, help="FASTA database file")
@click.option("--column", default="sequence", help="Column name for sequences (default: sequence)")
@click.option("--output", required=True, help="Output TSV file path")
def main(input, fasta, column, output):
    proteins = load_fasta(fasta)
    print(f"Loaded {len(proteins)} proteins from {fasta}")

    peptides = read_peptides_from_tsv(input, column=column)
    print(f"Read {len(peptides)} peptides from {input}")

    protein_peptides = map_peptides_to_proteins(peptides, proteins)
    print(f"Mapped peptides to {len(protein_peptides)} proteins")

    groups = build_protein_groups(protein_peptides, proteins)
    print(f"Built {len(groups)} protein groups")

    write_tsv(groups, output)
    print(f"Results written to {output}")


if __name__ == "__main__":
    main()
