"""
Library Coverage Estimator
===========================
Given a spectral library (TSV) and a FASTA proteome, estimate proteome
coverage: what fraction of theoretically digestible peptides appear in the
library, and what fraction of proteins have at least one peptide represented.

Uses pyopenms FASTAFile for reading the proteome, ProteaseDigestion for
in-silico digestion, and AASequence for sequence handling.

Usage
-----
    python library_coverage_estimator.py --library lib.tsv \
        --fasta proteome.fasta --enzyme Trypsin --output coverage.tsv
"""

import csv
from typing import Dict, List, Set, Tuple

import click
import pyopenms as oms


def read_library_peptides(library_path: str) -> Set[str]:
    """Read peptide sequences from a spectral library TSV.

    Expects a column named ``PeptideSequence`` or ``sequence``.

    Returns
    -------
    set
        Unique stripped peptide sequences.
    """
    peptides: Set[str] = set()
    with open(library_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            seq = row.get("PeptideSequence", row.get("sequence", "")).strip()
            if seq:
                # Strip modifications for matching
                aa = oms.AASequence.fromString(seq)
                peptides.add(aa.toUnmodifiedString())
    return peptides


def digest_fasta(
    fasta_path: str,
    enzyme: str = "Trypsin",
    missed_cleavages: int = 1,
    min_length: int = 6,
    max_length: int = 50,
) -> Tuple[Dict[str, List[str]], Set[str]]:
    """Digest a FASTA file and return per-protein peptides.

    Parameters
    ----------
    fasta_path:
        Path to FASTA file.
    enzyme:
        Enzyme name (default: Trypsin).
    missed_cleavages:
        Allowed missed cleavages (default: 1).
    min_length:
        Minimum peptide length.
    max_length:
        Maximum peptide length.

    Returns
    -------
    tuple
        (protein_peptides, all_peptides) where protein_peptides maps
        accession to list of peptide strings and all_peptides is the
        union of all digestible peptides.
    """
    entries: List[oms.FASTAEntry] = []
    fasta_file = oms.FASTAFile()
    fasta_file.load(fasta_path, entries)

    digester = oms.ProteaseDigestion()
    digester.setEnzyme(enzyme)
    digester.setMissedCleavages(missed_cleavages)

    protein_peptides: Dict[str, List[str]] = {}
    all_peptides: Set[str] = set()

    for entry in entries:
        accession = entry.identifier.split()[0] if entry.identifier else "unknown"
        aa_seq = oms.AASequence.fromString(entry.sequence)
        digest_result: List[oms.AASequence] = []
        digester.digest(aa_seq, digest_result, min_length, max_length)

        pep_strings = [p.toUnmodifiedString() for p in digest_result]
        protein_peptides[accession] = pep_strings
        all_peptides.update(pep_strings)

    return protein_peptides, all_peptides


def compute_coverage(
    library_peptides: Set[str],
    protein_peptides: Dict[str, List[str]],
    all_digestible: Set[str],
) -> Dict[str, object]:
    """Compute proteome coverage statistics.

    Parameters
    ----------
    library_peptides:
        Peptides present in the spectral library.
    protein_peptides:
        Per-protein peptide lists from digestion.
    all_digestible:
        Union of all digestible peptides.

    Returns
    -------
    dict
        Coverage metrics including peptide-level and protein-level fractions.
    """
    covered_peptides = library_peptides & all_digestible
    proteins_with_coverage = 0
    protein_details: List[Dict[str, object]] = []

    for acc, peps in protein_peptides.items():
        pep_set = set(peps)
        matched = pep_set & library_peptides
        has_coverage = len(matched) > 0
        if has_coverage:
            proteins_with_coverage += 1
        protein_details.append({
            "accession": acc,
            "total_peptides": len(pep_set),
            "library_peptides": len(matched),
            "coverage_fraction": len(matched) / len(pep_set) if pep_set else 0.0,
        })

    total_proteins = len(protein_peptides)
    return {
        "total_digestible_peptides": len(all_digestible),
        "library_peptides_in_proteome": len(covered_peptides),
        "peptide_coverage_fraction": (
            len(covered_peptides) / len(all_digestible) if all_digestible else 0.0
        ),
        "total_proteins": total_proteins,
        "proteins_with_library_peptide": proteins_with_coverage,
        "protein_coverage_fraction": (
            proteins_with_coverage / total_proteins if total_proteins else 0.0
        ),
        "protein_details": protein_details,
    }


@click.command(help="Estimate proteome coverage from a spectral library and FASTA.")
@click.option("--library", required=True, help="Spectral library TSV")
@click.option("--fasta", required=True, help="Proteome FASTA file")
@click.option("--enzyme", default="Trypsin", help="Enzyme name (default: Trypsin)")
@click.option("--missed-cleavages", type=int, default=1, help="Missed cleavages (default: 1)")
@click.option("--output", required=True, help="Output coverage TSV")
def main(library, fasta, enzyme, missed_cleavages, output) -> None:
    library_peps = read_library_peptides(library)
    protein_peps, all_peps = digest_fasta(fasta, enzyme, missed_cleavages)
    result = compute_coverage(library_peps, protein_peps, all_peps)

    with open(output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        # Summary
        writer.writerow(["metric", "value"])
        writer.writerow(["total_digestible_peptides", result["total_digestible_peptides"]])
        writer.writerow(["library_peptides_in_proteome", result["library_peptides_in_proteome"]])
        writer.writerow(["peptide_coverage_fraction", f"{result['peptide_coverage_fraction']:.4f}"])
        writer.writerow(["total_proteins", result["total_proteins"]])
        writer.writerow(["proteins_with_library_peptide", result["proteins_with_library_peptide"]])
        writer.writerow(["protein_coverage_fraction", f"{result['protein_coverage_fraction']:.4f}"])
        writer.writerow([])
        # Per-protein details
        writer.writerow(["accession", "total_peptides", "library_peptides", "coverage_fraction"])
        for det in result["protein_details"]:
            writer.writerow([
                det["accession"], det["total_peptides"],
                det["library_peptides"], f"{det['coverage_fraction']:.4f}",
            ])

    print(
        f"Peptide coverage: {result['library_peptides_in_proteome']}"
        f"/{result['total_digestible_peptides']}"
        f" ({result['peptide_coverage_fraction']:.1%})"
    )
    print(
        f"Protein coverage: {result['proteins_with_library_peptide']}"
        f"/{result['total_proteins']}"
        f" ({result['protein_coverage_fraction']:.1%})"
    )


if __name__ == "__main__":
    main()
