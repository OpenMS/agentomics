"""
Metapeptide LCA Assigner
=========================
Compute lowest common ancestor (LCA) taxonomy from peptide-protein mappings.

For each peptide, find all proteins it maps to, look up their taxonomy lineages,
and compute the LCA as the longest common prefix of all lineage strings.

Usage
-----
    python metapeptide_lca_assigner.py --peptides peptides.tsv --fasta metadb.fasta \
        --taxonomy lineage.tsv --output taxonomy.tsv
"""

import csv
import re
from typing import Dict, List, Set

import click
import pyopenms as oms


def load_fasta(fasta_path: str) -> Dict[str, str]:
    """Load a FASTA file into a dictionary mapping accession to sequence.

    Parameters
    ----------
    fasta_path:
        Path to the FASTA file.

    Returns
    -------
    dict
        Mapping of protein accession to amino acid sequence.
    """
    entries = []
    fasta_file = oms.FASTAFile()
    fasta_file.load(fasta_path, entries)

    proteins = {}
    for entry in entries:
        acc = entry.identifier.split()[0] if entry.identifier else ""
        proteins[acc] = entry.sequence
    return proteins


def load_taxonomy(taxonomy_path: str) -> Dict[str, List[str]]:
    """Load taxonomy lineage file.

    Parameters
    ----------
    taxonomy_path:
        Path to TSV file with columns: protein, lineage.
        Lineage is semicolon-separated taxonomy levels, e.g.
        "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia"

    Returns
    -------
    dict
        Mapping of protein accession to lineage list.
    """
    taxonomy: Dict[str, List[str]] = {}
    with open(taxonomy_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            acc = row["protein"]
            lineage_str = row["lineage"]
            lineage = [level.strip() for level in lineage_str.split(";") if level.strip()]
            taxonomy[acc] = lineage
    return taxonomy


def strip_modifications(sequence: str) -> str:
    """Remove modification annotations from a peptide sequence.

    Parameters
    ----------
    sequence:
        Peptide sequence with possible modifications.

    Returns
    -------
    str
        Clean amino acid sequence.
    """
    clean = re.sub(r"\[.*?\]", "", sequence)
    clean = re.sub(r"\(.*?\)", "", clean)
    return clean


def find_peptide_proteins(peptide: str, proteins: Dict[str, str]) -> Set[str]:
    """Find all proteins containing a given peptide sequence.

    Parameters
    ----------
    peptide:
        Peptide sequence.
    proteins:
        Protein accession to sequence mapping.

    Returns
    -------
    set
        Set of matching protein accessions.
    """
    clean = strip_modifications(peptide)
    matching = set()
    for acc, seq in proteins.items():
        if clean in seq:
            matching.add(acc)
    return matching


def compute_lca(lineages: List[List[str]]) -> List[str]:
    """Compute lowest common ancestor as the longest common prefix of lineages.

    Parameters
    ----------
    lineages:
        List of taxonomy lineage lists.

    Returns
    -------
    list
        LCA lineage (common prefix).
    """
    if not lineages:
        return []
    if len(lineages) == 1:
        return lineages[0]

    lca = []
    min_len = min(len(lin) for lin in lineages)
    for i in range(min_len):
        levels = {lin[i] for lin in lineages}
        if len(levels) == 1:
            lca.append(lineages[0][i])
        else:
            break
    return lca


def assign_lca_for_peptide(
    peptide: str,
    proteins: Dict[str, str],
    taxonomy: Dict[str, List[str]],
) -> Dict[str, object]:
    """Assign LCA taxonomy for a single peptide.

    Parameters
    ----------
    peptide:
        Peptide sequence.
    proteins:
        Protein accession to sequence mapping.
    taxonomy:
        Protein to lineage mapping.

    Returns
    -------
    dict
        Result with peptide, matched proteins, lineages, and LCA.
    """
    matched_prots = find_peptide_proteins(peptide, proteins)

    lineages = []
    for prot in sorted(matched_prots):
        if prot in taxonomy:
            lineages.append(taxonomy[prot])

    lca = compute_lca(lineages)

    return {
        "peptide": peptide,
        "matched_proteins": ";".join(sorted(matched_prots)) if matched_prots else "NONE",
        "num_proteins": len(matched_prots),
        "num_lineages": len(lineages),
        "lca": ";".join(lca) if lca else "unassigned",
        "lca_depth": len(lca),
        "taxonomic_specificity": lca[-1] if lca else "unassigned",
    }


def assign_lca_batch(
    peptides: List[str],
    proteins: Dict[str, str],
    taxonomy: Dict[str, List[str]],
) -> List[Dict[str, object]]:
    """Assign LCA for a batch of peptides.

    Parameters
    ----------
    peptides:
        List of peptide sequences.
    proteins:
        Protein accession to sequence mapping.
    taxonomy:
        Protein to lineage mapping.

    Returns
    -------
    list
        List of LCA assignment results.
    """
    results = []
    for pep in peptides:
        results.append(assign_lca_for_peptide(pep, proteins, taxonomy))
    return results


def read_peptides(peptides_path: str) -> List[str]:
    """Read peptides from TSV file with 'peptide' column."""
    peptides = []
    with open(peptides_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            peptides.append(row["peptide"])
    return peptides


def write_output(output_path: str, results: List[Dict[str, object]]) -> None:
    """Write LCA results to TSV."""
    if not results:
        return
    fieldnames = list(results[0].keys())
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


@click.command(help="Compute lowest common ancestor taxonomy from peptide-protein mappings.")
@click.option("--peptides", required=True, help="Peptides TSV file (peptide column)")
@click.option("--fasta", required=True, help="Meta-proteomics database FASTA file")
@click.option("--taxonomy", required=True, help="Taxonomy lineage TSV (protein, lineage)")
@click.option("--output", required=True, help="Output taxonomy TSV file")
def main(peptides, fasta, taxonomy, output):
    proteins = load_fasta(fasta)
    taxonomy_data = load_taxonomy(taxonomy)
    peptides_list = read_peptides(peptides)
    results = assign_lca_batch(peptides_list, proteins, taxonomy_data)
    write_output(output, results)

    n_assigned = sum(1 for r in results if r["lca"] != "unassigned")
    print(f"Total peptides: {len(results)}")
    print(f"Assigned LCA: {n_assigned}")
    print(f"Unassigned: {len(results) - n_assigned}")


if __name__ == "__main__":
    main()
