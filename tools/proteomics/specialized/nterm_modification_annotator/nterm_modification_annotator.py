"""
N-terminal Modification Annotator
===================================
Classify N-terminal peptides as protein N-terminus, signal peptide cleavage,
neo-N-terminus (internal cleavage), or unknown.

The tool maps each peptide to its protein, determines where it starts in the
protein sequence, and classifies the N-terminal type accordingly.

Usage
-----
    python nterm_modification_annotator.py --input nterm_peptides.tsv --fasta reference.fasta --output annotated.tsv
"""

import argparse
import csv
import re
import sys
from typing import Dict, List

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


# Typical signal peptide cleavage sites - residues that commonly precede signal peptide cleavage
SIGNAL_PEPTIDE_MAX_POS = 70  # Signal peptides are typically 15-70 residues


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


def get_clean_sequence(sequence: str) -> str:
    """Get clean amino acid sequence using AASequence.

    Parameters
    ----------
    sequence:
        Peptide sequence, possibly modified.

    Returns
    -------
    str
        Clean unmodified sequence.
    """
    try:
        aa = oms.AASequence.fromString(sequence)
        return aa.toUnmodifiedString()
    except Exception:
        clean = re.sub(r"\[.*?\]", "", sequence)
        clean = re.sub(r"\(.*?\)", "", clean)
        return clean


def detect_nterm_modification(sequence: str) -> str:
    """Detect N-terminal modifications from the peptide sequence string.

    Parameters
    ----------
    sequence:
        Original peptide sequence with modifications.

    Returns
    -------
    str
        Detected N-terminal modification or 'none'.
    """
    nterm_mods = {
        "Acetyl": "acetylation",
        "acetyl": "acetylation",
        ".(Acetyl)": "acetylation",
        "Formyl": "formylation",
        "formyl": "formylation",
        "Carbamyl": "carbamylation",
        "TMT": "TMT-label",
        "iTRAQ": "iTRAQ-label",
        "Dimethyl": "dimethylation",
    }

    for mod_key, mod_name in nterm_mods.items():
        if mod_key in sequence:
            return mod_name
    return "none"


def find_peptide_start(peptide_clean: str, protein_seq: str) -> int:
    """Find the 0-based start position of a peptide in a protein sequence.

    Parameters
    ----------
    peptide_clean:
        Clean peptide sequence.
    protein_seq:
        Full protein sequence.

    Returns
    -------
    int
        0-based start position, or -1 if not found.
    """
    return protein_seq.find(peptide_clean)


def classify_nterm_type(
    start_pos: int,
    protein_seq: str,
    signal_peptide_sites: Dict[str, int] = None,
    protein_id: str = "",
) -> str:
    """Classify the N-terminal peptide type.

    Parameters
    ----------
    start_pos:
        0-based start position of peptide in protein.
    protein_seq:
        Full protein sequence.
    signal_peptide_sites:
        Optional dict of protein_id -> signal peptide cleavage position (0-based).
    protein_id:
        Protein accession for signal peptide lookup.

    Returns
    -------
    str
        Classification: 'protein_nterm', 'met_removal', 'signal_peptide', 'neo_nterm', or 'unmapped'.
    """
    if start_pos < 0:
        return "unmapped"

    if start_pos == 0:
        return "protein_nterm"

    # Methionine removal: peptide starts at position 1 and protein starts with M
    if start_pos == 1 and protein_seq[0] == "M":
        return "met_removal"

    # Signal peptide cleavage: check known sites or heuristic
    if signal_peptide_sites and protein_id in signal_peptide_sites:
        sp_site = signal_peptide_sites[protein_id]
        if start_pos == sp_site:
            return "signal_peptide"

    # Heuristic: if start position is within typical signal peptide range
    # and preceded by small/neutral residues (A, G, S are common at -1 position)
    if 15 <= start_pos <= SIGNAL_PEPTIDE_MAX_POS:
        preceding_aa = protein_seq[start_pos - 1]
        if preceding_aa in "AGS":
            return "signal_peptide_candidate"

    return "neo_nterm"


def annotate_nterm_peptides(
    rows: List[Dict[str, str]],
    proteins: Dict[str, str],
    signal_peptide_sites: Dict[str, int] = None,
) -> List[Dict[str, str]]:
    """Annotate a batch of N-terminal peptides.

    Parameters
    ----------
    rows:
        List of dicts with keys: peptide, protein.
    proteins:
        Protein accession to sequence mapping.
    signal_peptide_sites:
        Optional dict of known signal peptide cleavage sites.

    Returns
    -------
    list
        Annotated rows with added classification columns.
    """
    results = []

    for row in rows:
        peptide_raw = row["peptide"]
        protein_id = row["protein"]
        clean_seq = get_clean_sequence(peptide_raw)

        new_row = dict(row)
        new_row["clean_sequence"] = clean_seq
        new_row["nterm_modification"] = detect_nterm_modification(peptide_raw)

        if protein_id not in proteins:
            new_row["start_position"] = -1
            new_row["nterm_type"] = "unmapped"
            results.append(new_row)
            continue

        protein_seq = proteins[protein_id]
        start_pos = find_peptide_start(clean_seq, protein_seq)
        nterm_type = classify_nterm_type(start_pos, protein_seq, signal_peptide_sites, protein_id)

        new_row["start_position"] = start_pos
        new_row["nterm_type"] = nterm_type
        results.append(new_row)

    return results


def compute_summary(results: List[Dict[str, str]]) -> Dict[str, int]:
    """Compute summary counts by N-terminal type.

    Parameters
    ----------
    results:
        List of annotated result dicts.

    Returns
    -------
    dict
        Counts per nterm_type.
    """
    summary: Dict[str, int] = {}
    for r in results:
        nt = r["nterm_type"]
        summary[nt] = summary.get(nt, 0) + 1
    return summary


def read_input(input_path: str) -> List[Dict[str, str]]:
    """Read N-terminal peptides TSV file."""
    rows = []
    with open(input_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def write_output(output_path: str, results: List[Dict[str, str]]) -> None:
    """Write annotated results to TSV."""
    if not results:
        return
    fieldnames = list(results[0].keys())
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


def main():
    parser = argparse.ArgumentParser(
        description="Classify N-terminal peptides as protein N-term, signal peptide, neo-N-term, etc."
    )
    parser.add_argument("--input", required=True, help="Input N-terminal peptides TSV file")
    parser.add_argument("--fasta", required=True, help="Reference proteome FASTA file")
    parser.add_argument("--output", required=True, help="Output annotated TSV file")
    args = parser.parse_args()

    proteins = load_fasta(args.fasta)
    rows = read_input(args.input)
    results = annotate_nterm_peptides(rows, proteins)
    write_output(args.output, results)

    summary = compute_summary(results)
    print(f"Total peptides: {len(results)}")
    for ntype, count in sorted(summary.items()):
        print(f"  {ntype}: {count}")


if __name__ == "__main__":
    main()
