"""
Crosslink Classifier
====================
Classify crosslinks as intra-protein, inter-protein, or monolink based on
peptide-to-protein mappings from a FASTA database.

Usage
-----
    python xl_link_classifier.py --crosslinks links.tsv --fasta proteome.fasta --output classified.tsv
"""

import argparse
import csv
import sys
from typing import Dict, List, Set

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


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


def find_peptide_proteins(peptide: str, proteins: Dict[str, str]) -> Set[str]:
    """Find all proteins that contain a given peptide sequence.

    Parameters
    ----------
    peptide:
        Peptide amino acid sequence (unmodified).
    proteins:
        Protein accession to sequence mapping.

    Returns
    -------
    set
        Set of protein accessions containing the peptide.
    """
    matching = set()
    # Strip modifications for matching
    clean = strip_modifications(peptide)
    for acc, seq in proteins.items():
        if clean in seq:
            matching.add(acc)
    return matching


def strip_modifications(sequence: str) -> str:
    """Remove modification annotations from a peptide sequence.

    Parameters
    ----------
    sequence:
        Peptide sequence possibly containing bracket/parenthesis modifications.

    Returns
    -------
    str
        Clean amino acid sequence.
    """
    import re
    clean = re.sub(r"\[.*?\]", "", sequence)
    clean = re.sub(r"\(.*?\)", "", clean)
    return clean


def classify_crosslink(
    peptide1: str,
    peptide2: str,
    proteins: Dict[str, str],
) -> Dict[str, object]:
    """Classify a crosslink as intra-protein, inter-protein, or monolink.

    Parameters
    ----------
    peptide1:
        First peptide sequence.
    peptide2:
        Second peptide sequence (empty string for monolinks).
    proteins:
        Protein accession to sequence mapping.

    Returns
    -------
    dict
        Classification result with link_type, proteins1, proteins2.
    """
    # Monolink: second peptide is empty or missing
    if not peptide2 or peptide2.strip() == "" or peptide2.strip() == "-":
        prots1 = find_peptide_proteins(peptide1, proteins)
        return {
            "peptide1": peptide1,
            "peptide2": peptide2,
            "link_type": "monolink",
            "proteins1": ";".join(sorted(prots1)) if prots1 else "UNMAPPED",
            "proteins2": "",
            "shared_proteins": "",
        }

    prots1 = find_peptide_proteins(peptide1, proteins)
    prots2 = find_peptide_proteins(peptide2, proteins)

    shared = prots1 & prots2

    if not prots1 or not prots2:
        link_type = "unknown"
    elif shared:
        # Both peptides map to at least one common protein
        if prots1 == prots2 and len(prots1) == 1:
            link_type = "intra-protein"
        elif shared:
            # Could be intra on shared protein(s), but also maps to others
            link_type = "intra-protein"
        else:
            link_type = "inter-protein"
    else:
        link_type = "inter-protein"

    return {
        "peptide1": peptide1,
        "peptide2": peptide2,
        "link_type": link_type,
        "proteins1": ";".join(sorted(prots1)) if prots1 else "UNMAPPED",
        "proteins2": ";".join(sorted(prots2)) if prots2 else "UNMAPPED",
        "shared_proteins": ";".join(sorted(shared)) if shared else "",
    }


def classify_crosslinks(
    crosslinks: List[Dict[str, str]],
    proteins: Dict[str, str],
) -> List[Dict[str, object]]:
    """Classify a batch of crosslinks.

    Parameters
    ----------
    crosslinks:
        List of dicts with keys: peptide1, peptide2.
    proteins:
        Protein accession to sequence mapping.

    Returns
    -------
    list
        List of classification result dicts.
    """
    results = []
    for xl in crosslinks:
        pep1 = xl.get("peptide1", "")
        pep2 = xl.get("peptide2", "")
        result = classify_crosslink(pep1, pep2, proteins)
        # Preserve extra columns from input
        for key, val in xl.items():
            if key not in result:
                result[key] = val
        results.append(result)
    return results


def compute_summary(results: List[Dict[str, object]]) -> Dict[str, int]:
    """Compute summary counts by link type.

    Parameters
    ----------
    results:
        List of classification result dicts.

    Returns
    -------
    dict
        Counts per link type.
    """
    summary: Dict[str, int] = {}
    for r in results:
        lt = str(r["link_type"])
        summary[lt] = summary.get(lt, 0) + 1
    return summary


def read_crosslinks(crosslinks_path: str) -> List[Dict[str, str]]:
    """Read crosslinks TSV file."""
    rows = []
    with open(crosslinks_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def write_output(output_path: str, results: List[Dict[str, object]]) -> None:
    """Write classified crosslinks to TSV."""
    if not results:
        return
    fieldnames = list(results[0].keys())
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


def main():
    parser = argparse.ArgumentParser(
        description="Classify crosslinks as intra/inter-protein or monolink."
    )
    parser.add_argument("--crosslinks", required=True, help="Crosslinks TSV file")
    parser.add_argument("--fasta", required=True, help="Proteome FASTA file")
    parser.add_argument("--output", required=True, help="Output classified TSV file")
    args = parser.parse_args()

    proteins = load_fasta(args.fasta)
    crosslinks = read_crosslinks(args.crosslinks)
    results = classify_crosslinks(crosslinks, proteins)
    write_output(args.output, results)

    summary = compute_summary(results)
    print(f"Total crosslinks: {len(results)}")
    for lt, count in sorted(summary.items()):
        print(f"  {lt}: {count}")


if __name__ == "__main__":
    main()
