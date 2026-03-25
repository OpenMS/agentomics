"""
Semi-Tryptic Peptide Finder
============================
Classify peptides as fully tryptic, semi-tryptic, or non-tryptic by checking
whether their N- and C-terminal cleavage sites match the enzyme specificity.

Usage
-----
    python semi_tryptic_peptide_finder.py --input peptides.tsv --fasta db.fasta --enzyme Trypsin --output classified.tsv
"""

import argparse
import csv
import sys
from typing import List

try:
    import pyopenms as oms
except ImportError:
    sys.exit(
        "pyopenms is required. Install it with:  pip install pyopenms"
    )


def load_fasta(fasta_path: str) -> dict:
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


def digest_protein(protein_sequence: str, enzyme: str = "Trypsin", missed_cleavages: int = 2) -> List[str]:
    """Digest a protein sequence with the given enzyme.

    Parameters
    ----------
    protein_sequence:
        Full protein amino acid sequence.
    enzyme:
        Enzyme name (default: Trypsin).
    missed_cleavages:
        Number of allowed missed cleavages.

    Returns
    -------
    list
        List of fully tryptic peptide strings.
    """
    digestion = oms.ProteaseDigestion()
    digestion.setEnzyme(enzyme)
    digestion.setMissedCleavages(missed_cleavages)

    aa_seq = oms.AASequence.fromString(protein_sequence)
    result = []
    digestion.digest(aa_seq, result)
    return [str(pep) for pep in result]


def classify_peptide(peptide: str, protein_sequence: str, enzyme: str = "Trypsin") -> str:
    """Classify a peptide as fully_tryptic, semi_tryptic, or non_tryptic.

    Parameters
    ----------
    peptide:
        Peptide sequence to classify.
    protein_sequence:
        Parent protein sequence.
    enzyme:
        Enzyme name.

    Returns
    -------
    str
        Classification: 'fully_tryptic', 'semi_tryptic', or 'non_tryptic'.
    """
    # Find peptide in protein
    pos = protein_sequence.find(peptide)
    if pos == -1:
        return "not_found"

    # Check N-terminal cleavage (trypsin: after K or R, or protein N-term)
    n_term_ok = False
    if pos == 0:
        n_term_ok = True
    elif enzyme == "Trypsin" and protein_sequence[pos - 1] in ("K", "R"):
        # Check for proline rule: trypsin does not cleave before P
        if peptide[0] != "P":
            n_term_ok = True

    # Check C-terminal cleavage (trypsin: peptide ends with K or R, or protein C-term)
    c_term_ok = False
    end_pos = pos + len(peptide)
    if end_pos == len(protein_sequence):
        c_term_ok = True
    elif enzyme == "Trypsin" and peptide[-1] in ("K", "R"):
        c_term_ok = True

    if n_term_ok and c_term_ok:
        return "fully_tryptic"
    elif n_term_ok or c_term_ok:
        return "semi_tryptic"
    else:
        return "non_tryptic"


def classify_peptides_against_fasta(
    peptides: List[str],
    proteins: dict,
    enzyme: str = "Trypsin",
) -> List[dict]:
    """Classify a list of peptides against protein sequences.

    Parameters
    ----------
    peptides:
        List of peptide sequences.
    proteins:
        Dict mapping accession to protein sequence.
    enzyme:
        Enzyme name.

    Returns
    -------
    list
        List of dicts with sequence, classification, and matched protein.
    """
    results = []
    for pep in peptides:
        pep = pep.strip()
        if not pep:
            continue
        best_class = "not_found"
        matched_protein = ""
        for acc, prot_seq in proteins.items():
            classification = classify_peptide(pep, prot_seq, enzyme)
            if classification == "fully_tryptic":
                best_class = classification
                matched_protein = acc
                break
            elif classification == "semi_tryptic" and best_class != "fully_tryptic":
                best_class = classification
                matched_protein = acc
            elif classification == "non_tryptic" and best_class == "not_found":
                best_class = classification
                matched_protein = acc

        aa_seq = oms.AASequence.fromString(pep)
        results.append({
            "sequence": pep,
            "length": aa_seq.size(),
            "classification": best_class,
            "protein": matched_protein,
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
    """Write classification results to TSV.

    Parameters
    ----------
    results:
        List of result dicts.
    output_path:
        Output file path.
    """
    fieldnames = ["sequence", "length", "classification", "protein"]
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in results:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(
        description="Classify peptides as fully/semi/non-tryptic."
    )
    parser.add_argument("--input", required=True, help="Input TSV with peptide sequences")
    parser.add_argument("--fasta", required=True, help="FASTA database file")
    parser.add_argument("--enzyme", default="Trypsin", help="Enzyme name (default: Trypsin)")
    parser.add_argument("--column", default="sequence", help="Column name for sequences (default: sequence)")
    parser.add_argument("--output", required=True, help="Output TSV file path")
    args = parser.parse_args()

    proteins = load_fasta(args.fasta)
    print(f"Loaded {len(proteins)} proteins from {args.fasta}")

    peptides = read_peptides_from_tsv(args.input, column=args.column)
    print(f"Read {len(peptides)} peptides from {args.input}")

    results = classify_peptides_against_fasta(peptides, proteins, enzyme=args.enzyme)

    counts = {}
    for r in results:
        counts[r["classification"]] = counts.get(r["classification"], 0) + 1
    for cls, cnt in sorted(counts.items()):
        print(f"  {cls}: {cnt}")

    write_tsv(results, args.output)
    print(f"Results written to {args.output}")


if __name__ == "__main__":
    main()
