"""
Cleavage Site Profiler
======================
From neo-N-terminal peptides, extract P4-P4' windows around cleavage sites and compute
position-specific amino acid frequencies.

The P4-P4' nomenclature (Schechter & Berger) describes residues surrounding a cleavage
site: P4-P3-P2-P1 | P1'-P2'-P3'-P4' where | is the cleavage site.

Usage
-----
    python cleavage_site_profiler.py --input neo_nterm.tsv --fasta reference.fasta --window 4 --output profile.tsv
"""

import csv
from collections import Counter
from typing import Dict, List, Tuple

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
        import re
        clean = re.sub(r"\[.*?\]", "", sequence)
        clean = re.sub(r"\(.*?\)", "", clean)
        return clean


def find_cleavage_position(peptide: str, protein: str, protein_seq: str) -> int:
    """Find the cleavage site position in the protein sequence.

    The cleavage site is the position just before the neo-N-terminal peptide starts.

    Parameters
    ----------
    peptide:
        Neo-N-terminal peptide sequence (clean).
    protein:
        Protein accession (unused, kept for API consistency).
    protein_seq:
        Full protein sequence.

    Returns
    -------
    int
        0-based position of the cleavage site (P1 position), or -1 if not found.
    """
    idx = protein_seq.find(peptide)
    if idx <= 0:
        return -1
    # Cleavage site is between position idx-1 (P1) and idx (P1')
    return idx


def extract_cleavage_window(protein_seq: str, cleavage_pos: int, window: int = 4) -> str:
    """Extract P(window)...P1-P1'...P(window)' window around a cleavage site.

    Parameters
    ----------
    protein_seq:
        Full protein sequence.
    cleavage_pos:
        0-based position of P1' (first residue of neo-N-terminal peptide).
    window:
        Number of residues on each side of the cleavage (default 4 for P4-P4').

    Returns
    -------
    str
        Window of length 2*window, padded with '_' at termini.
    """
    result = []
    # P-side: positions cleavage_pos-window to cleavage_pos-1 (P4..P1)
    for i in range(cleavage_pos - window, cleavage_pos):
        if i < 0 or i >= len(protein_seq):
            result.append("_")
        else:
            result.append(protein_seq[i])
    # P'-side: positions cleavage_pos to cleavage_pos+window-1 (P1'..P4')
    for i in range(cleavage_pos, cleavage_pos + window):
        if i < 0 or i >= len(protein_seq):
            result.append("_")
        else:
            result.append(protein_seq[i])
    return "".join(result)


def compute_position_frequencies(
    windows: List[str], window: int = 4
) -> Dict[str, Dict[str, float]]:
    """Compute position-specific amino acid frequencies.

    Parameters
    ----------
    windows:
        List of cleavage site window strings.
    window:
        Window size on each side.

    Returns
    -------
    dict
        Mapping of position label (P4..P1, P1'..P4') to amino acid frequency dict.
    """
    total = len(windows)
    if total == 0:
        return {}

    labels = []
    for i in range(window, 0, -1):
        labels.append(f"P{i}")
    for i in range(1, window + 1):
        labels.append(f"P{i}'")

    frequencies: Dict[str, Dict[str, float]] = {}
    for pos_idx, label in enumerate(labels):
        counter: Counter = Counter()
        for w in windows:
            if pos_idx < len(w):
                counter[w[pos_idx]] += 1
        frequencies[label] = {aa: count / total for aa, count in counter.most_common()}

    return frequencies


def process_neo_nterm_peptides(
    rows: List[Dict[str, str]],
    proteins: Dict[str, str],
    window: int = 4,
) -> Tuple[List[Dict[str, str]], List[str]]:
    """Process neo-N-terminal peptides to extract cleavage windows.

    Parameters
    ----------
    rows:
        List of dicts with keys: peptide, protein.
    proteins:
        Protein accession to sequence mapping.
    window:
        Cleavage window size.

    Returns
    -------
    tuple
        (result_rows with added cleavage_window, list of valid windows)
    """
    results = []
    valid_windows = []

    for row in rows:
        peptide_raw = row["peptide"]
        protein_id = row["protein"]
        clean_seq = get_clean_sequence(peptide_raw)

        new_row = dict(row)
        new_row["clean_sequence"] = clean_seq

        if protein_id not in proteins:
            new_row["cleavage_window"] = "_" * (2 * window)
            new_row["cleavage_found"] = "NO"
            results.append(new_row)
            continue

        protein_seq = proteins[protein_id]
        cleavage_pos = find_cleavage_position(clean_seq, protein_id, protein_seq)

        if cleavage_pos < 0:
            new_row["cleavage_window"] = "_" * (2 * window)
            new_row["cleavage_found"] = "NO"
        else:
            w = extract_cleavage_window(protein_seq, cleavage_pos, window)
            new_row["cleavage_window"] = w
            new_row["cleavage_found"] = "YES"
            valid_windows.append(w)

        results.append(new_row)

    return results, valid_windows


def read_input(input_path: str) -> List[Dict[str, str]]:
    """Read neo-N-terminal peptides TSV file."""
    rows = []
    with open(input_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def write_output(
    output_path: str,
    result_rows: List[Dict[str, str]],
    frequencies: Dict[str, Dict[str, float]],
) -> None:
    """Write cleavage profiles to output files."""
    if result_rows:
        fieldnames = list(result_rows[0].keys())
        with open(output_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(result_rows)

    freq_path = output_path.replace(".tsv", "_frequencies.tsv")
    with open(freq_path, "w", newline="") as f:
        f.write("position\tamino_acid\tfrequency\n")
        for pos_label in frequencies:
            for aa, freq in sorted(frequencies[pos_label].items(), key=lambda x: -x[1]):
                f.write(f"{pos_label}\t{aa}\t{freq:.4f}\n")


@click.command(help="Profile cleavage sites from neo-N-terminal peptides.")
@click.option("--input", "input", required=True, help="Neo-N-terminal peptides TSV file")
@click.option("--fasta", required=True, help="Reference proteome FASTA file")
@click.option("--window", type=int, default=4, help="Window size on each side (default: 4)")
@click.option("--output", required=True, help="Output profile TSV file")
def main(input, fasta, window, output):
    proteins = load_fasta(fasta)
    rows = read_input(input)
    result_rows, valid_windows = process_neo_nterm_peptides(rows, proteins, window)
    frequencies = compute_position_frequencies(valid_windows, window)
    write_output(output, result_rows, frequencies)

    print(f"Total peptides: {len(result_rows)}")
    print(f"Cleavage sites found: {len(valid_windows)}")
    print(f"Window size: P{window}-P{window}'")


if __name__ == "__main__":
    main()
