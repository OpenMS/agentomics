"""
Phospho Motif Analyzer
======================
Extract +/-7 amino acid windows around phosphosites, compute position-specific frequencies.

Given a list of phosphosites (peptide, protein, site position) and a FASTA proteome,
this tool extracts the surrounding sequence window and computes amino acid frequency
at each position relative to the phosphosite.

Usage
-----
    python phospho_motif_analyzer.py --input phosphosites.tsv --fasta proteome.fasta --window 7 --output motifs.tsv
"""

import csv
from collections import Counter
from typing import Dict, List, Tuple

import click
import pyopenms as oms


def load_fasta(fasta_path: str) -> Dict[str, str]:
    """Load a FASTA file into a dictionary mapping protein accession to sequence.

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
        # Use the identifier (first word of description or accession)
        acc = entry.identifier.split()[0] if entry.identifier else ""
        proteins[acc] = entry.sequence
    return proteins


def extract_window(protein_seq: str, site_pos: int, window: int = 7) -> str:
    """Extract a window of amino acids around a site position.

    Parameters
    ----------
    protein_seq:
        Full protein sequence.
    site_pos:
        1-based position of the phosphosite in the protein.
    window:
        Number of residues on each side (default 7).

    Returns
    -------
    str
        Window sequence of length 2*window+1, padded with '_' if near terminus.
    """
    idx = site_pos - 1  # Convert to 0-based
    start = idx - window
    end = idx + window + 1

    result = []
    for i in range(start, end):
        if i < 0 or i >= len(protein_seq):
            result.append("_")
        else:
            result.append(protein_seq[i])
    return "".join(result)


def validate_peptide(sequence: str) -> bool:
    """Validate a peptide sequence using AASequence.

    Parameters
    ----------
    sequence:
        Peptide sequence string.

    Returns
    -------
    bool
        True if parseable.
    """
    try:
        oms.AASequence.fromString(sequence)
        return True
    except Exception:
        return False


def extract_motif_windows(
    rows: List[Dict[str, str]],
    proteins: Dict[str, str],
    window: int = 7,
) -> List[Dict[str, str]]:
    """Extract motif windows for all phosphosite rows.

    Parameters
    ----------
    rows:
        List of dicts with keys: peptide, protein, site (1-based position).
    proteins:
        Protein accession to sequence mapping.
    window:
        Window size on each side.

    Returns
    -------
    list
        List of dicts with added 'motif_window' key.
    """
    results = []
    for row in rows:
        protein_id = row["protein"]
        site_pos = int(row["site"])
        new_row = dict(row)

        if protein_id in proteins:
            motif = extract_window(proteins[protein_id], site_pos, window)
            new_row["motif_window"] = motif
        else:
            new_row["motif_window"] = "_" * (2 * window + 1)

        new_row["valid_peptide"] = str(validate_peptide(row["peptide"]))
        results.append(new_row)

    return results


def compute_position_frequencies(
    windows: List[str], window: int = 7
) -> Dict[int, Dict[str, float]]:
    """Compute position-specific amino acid frequencies from a set of motif windows.

    Parameters
    ----------
    windows:
        List of motif window strings (all same length).
    window:
        Window size used to extract motifs.

    Returns
    -------
    dict
        Mapping of relative position (-window to +window) to amino acid frequency dict.
    """
    total = len(windows)
    if total == 0:
        return {}

    width = 2 * window + 1
    frequencies: Dict[int, Dict[str, float]] = {}

    for pos in range(width):
        rel_pos = pos - window
        counter: Counter = Counter()
        for w in windows:
            if pos < len(w):
                counter[w[pos]] += 1
        frequencies[rel_pos] = {aa: count / total for aa, count in counter.most_common()}

    return frequencies


def format_frequencies(frequencies: Dict[int, Dict[str, float]]) -> List[Tuple[int, str, float]]:
    """Flatten frequency dict into a list of (position, amino_acid, frequency) tuples.

    Parameters
    ----------
    frequencies:
        Position-specific frequency dict.

    Returns
    -------
    list
        List of (position, amino_acid, frequency) tuples.
    """
    result = []
    for pos in sorted(frequencies.keys()):
        for aa, freq in sorted(frequencies[pos].items(), key=lambda x: -x[1]):
            result.append((pos, aa, freq))
    return result


def read_input(input_path: str) -> List[Dict[str, str]]:
    """Read phosphosite TSV input file."""
    rows = []
    with open(input_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def write_output(
    output_path: str,
    motif_rows: List[Dict[str, str]],
    frequencies: Dict[int, Dict[str, float]],
) -> None:
    """Write motif windows and frequencies to output files."""
    if motif_rows:
        fieldnames = list(motif_rows[0].keys())
        with open(output_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(motif_rows)

    freq_path = output_path.replace(".tsv", "_frequencies.tsv")
    flat = format_frequencies(frequencies)
    with open(freq_path, "w", newline="") as f:
        f.write("position\tamino_acid\tfrequency\n")
        for pos, aa, freq in flat:
            f.write(f"{pos}\t{aa}\t{freq:.4f}\n")


@click.command(help="Extract motif windows around phosphosites and compute position-specific frequencies.")
@click.option("--input", "input", required=True, help="Input phosphosites TSV file")
@click.option("--fasta", required=True, help="Proteome FASTA file")
@click.option("--window", type=int, default=7, help="Window size on each side (default: 7)")
@click.option("--output", required=True, help="Output motifs TSV file")
def main(input, fasta, window, output):
    proteins = load_fasta(fasta)
    rows = read_input(input)
    motif_rows = extract_motif_windows(rows, proteins, window)
    windows_list = [r["motif_window"] for r in motif_rows]
    frequencies = compute_position_frequencies(windows_list, window)
    write_output(output, motif_rows, frequencies)

    print(f"Processed {len(motif_rows)} phosphosites")
    print(f"Window size: +/-{window} residues")
    print(f"Output written to {output}")


if __name__ == "__main__":
    main()
