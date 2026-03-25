"""
Metapeptide Function Aggregator
================================
Aggregate GO/KEGG functional annotations from peptide-to-protein mappings.
Given identified peptides with their protein assignments and a separate
annotation file mapping proteins to functional terms, the tool aggregates
term counts and computes peptide-level functional profiles.

Useful for metaproteomics where functional characterization of the
community is derived from identified peptides.

Usage
-----
    python metapeptide_function_aggregator.py --peptides identified.tsv \
        --annotations go_terms.tsv --output function.tsv
"""

import argparse
import csv
import sys
from collections import Counter, defaultdict
from typing import Dict, List, Set, Tuple

try:
    import pyopenms as oms  # noqa: F401
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def load_peptide_protein_map(peptides_path: str) -> Dict[str, Set[str]]:
    """Load peptide-to-protein mappings from a TSV.

    Expects columns ``peptide`` and ``protein``.  A peptide can map to
    multiple proteins (one row per mapping or semicolon-separated proteins).

    Returns
    -------
    dict
        Mapping of peptide sequence to set of protein accessions.
    """
    pep_to_prot: Dict[str, Set[str]] = defaultdict(set)
    with open(peptides_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            peptide = row.get("peptide", "").strip()
            proteins_raw = row.get("protein", "").strip()
            if not peptide or not proteins_raw:
                continue
            for prot in proteins_raw.split(";"):
                prot = prot.strip()
                if prot:
                    pep_to_prot[peptide].add(prot)
    return dict(pep_to_prot)


def load_annotations(annotations_path: str) -> Dict[str, List[Tuple[str, str]]]:
    """Load protein-to-function annotation mappings.

    Expects columns ``protein``, ``term_id``, and ``term_name``.

    Returns
    -------
    dict
        Mapping of protein accession to list of (term_id, term_name) tuples.
    """
    prot_to_terms: Dict[str, List[Tuple[str, str]]] = defaultdict(list)
    with open(annotations_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            protein = row.get("protein", "").strip()
            term_id = row.get("term_id", "").strip()
            term_name = row.get("term_name", "").strip()
            if protein and term_id:
                prot_to_terms[protein].append((term_id, term_name))
    return dict(prot_to_terms)


def aggregate_functions(
    pep_to_prot: Dict[str, Set[str]],
    prot_to_terms: Dict[str, List[Tuple[str, str]]],
) -> Tuple[List[Dict[str, object]], Counter]:
    """Aggregate functional annotations from peptide-protein-term mappings.

    For each peptide, collect all functional terms from all mapped proteins.
    Also compute global term frequency counts.

    Parameters
    ----------
    pep_to_prot:
        Peptide-to-protein mappings.
    prot_to_terms:
        Protein-to-functional-term mappings.

    Returns
    -------
    tuple
        (peptide_annotations, term_counts) where peptide_annotations is a list
        of dicts with ``peptide``, ``proteins``, ``terms``, and term_counts is
        a Counter of term_id occurrences across all peptides.
    """
    peptide_annotations: List[Dict[str, object]] = []
    term_counts: Counter = Counter()

    for peptide, proteins in sorted(pep_to_prot.items()):
        terms_seen: Set[str] = set()
        term_details: List[Tuple[str, str]] = []

        for prot in proteins:
            if prot in prot_to_terms:
                for term_id, term_name in prot_to_terms[prot]:
                    if term_id not in terms_seen:
                        terms_seen.add(term_id)
                        term_details.append((term_id, term_name))
                        term_counts[term_id] += 1

        peptide_annotations.append({
            "peptide": peptide,
            "n_proteins": len(proteins),
            "proteins": ";".join(sorted(proteins)),
            "n_terms": len(term_details),
            "terms": ";".join(f"{tid}:{tname}" for tid, tname in term_details),
        })

    return peptide_annotations, term_counts


def summarize_terms(
    term_counts: Counter,
    prot_to_terms: Dict[str, List[Tuple[str, str]]],
) -> List[Dict[str, object]]:
    """Build a summary table of term frequencies.

    Returns
    -------
    list of dict
        Sorted by count descending, with ``term_id``, ``term_name``, ``count``.
    """
    # Build term_id -> term_name lookup
    id_to_name: Dict[str, str] = {}
    for terms in prot_to_terms.values():
        for tid, tname in terms:
            if tid not in id_to_name:
                id_to_name[tid] = tname

    summary: List[Dict[str, object]] = []
    for tid, count in term_counts.most_common():
        summary.append({
            "term_id": tid,
            "term_name": id_to_name.get(tid, ""),
            "peptide_count": count,
        })
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Aggregate GO/KEGG annotations from peptide-protein mappings."
    )
    parser.add_argument(
        "--peptides", required=True,
        help="TSV with 'peptide' and 'protein' columns",
    )
    parser.add_argument(
        "--annotations", required=True,
        help="TSV with 'protein', 'term_id', 'term_name' columns",
    )
    parser.add_argument("--output", required=True, help="Output functional aggregation TSV")
    args = parser.parse_args()

    pep_to_prot = load_peptide_protein_map(args.peptides)
    prot_to_terms = load_annotations(args.annotations)

    if not pep_to_prot:
        sys.exit("No peptide-protein mappings found.")

    pep_annots, term_counts = aggregate_functions(pep_to_prot, prot_to_terms)
    term_summary = summarize_terms(term_counts, prot_to_terms)

    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        # Peptide-level annotations
        writer.writerow(["peptide", "n_proteins", "proteins", "n_terms", "terms"])
        for pa in pep_annots:
            writer.writerow([
                pa["peptide"], pa["n_proteins"], pa["proteins"],
                pa["n_terms"], pa["terms"],
            ])
        writer.writerow([])
        # Term summary
        writer.writerow(["term_id", "term_name", "peptide_count"])
        for ts in term_summary:
            writer.writerow([ts["term_id"], ts["term_name"], ts["peptide_count"]])

    annotated = sum(1 for pa in pep_annots if pa["n_terms"] > 0)
    print(f"Annotated {annotated}/{len(pep_annots)} peptides with functional terms")
    print(f"Unique terms: {len(term_summary)} -> {args.output}")


if __name__ == "__main__":
    main()
