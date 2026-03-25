"""
Phospho Enrichment QC
=====================
Compute phospho-enrichment efficiency and pSer/pThr/pTyr ratios from search results.

Analyzes peptide search results to determine what fraction of identified peptides
carry phosphorylation modifications, and breaks down the phosphorylation by
residue type (Ser, Thr, Tyr).

Usage
-----
    python phospho_enrichment_qc.py --input search_results.tsv --output enrichment.tsv
"""

import csv
import re
import sys
from typing import Dict, List, Tuple

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def is_phosphopeptide(sequence: str) -> bool:
    """Check if a peptide sequence contains phosphorylation modifications.

    Parameters
    ----------
    sequence:
        Peptide sequence, possibly with modifications in bracket notation.

    Returns
    -------
    bool
        True if the sequence contains phosphorylation.
    """
    phospho_patterns = ["Phospho", "phospho", "(Phospho)", "[80]", "[79.966]"]
    for pat in phospho_patterns:
        if pat in sequence:
            return True
    return False


def count_phospho_residues(sequence: str) -> Dict[str, int]:
    """Count phosphorylated Ser, Thr, and Tyr residues in a modified peptide sequence.

    Parameters
    ----------
    sequence:
        Modified peptide sequence string.

    Returns
    -------
    dict
        Counts with keys 'pSer', 'pThr', 'pTyr'.
    """
    counts = {"pSer": 0, "pThr": 0, "pTyr": 0}

    # Pattern: residue followed by modification indicating phospho
    # Matches S(Phospho), T(Phospho), Y(Phospho) or S[80], T[80], Y[80] etc.
    phospho_mods = [r"S\(Phospho\)", r"T\(Phospho\)", r"Y\(Phospho\)",
                    r"S\[80\]", r"T\[80\]", r"Y\[80\]",
                    r"S\[79\.966\]", r"T\[79\.966\]", r"Y\[79\.966\]"]

    for pattern in phospho_mods[:3]:
        counts["pSer"] += len(re.findall(pattern, sequence))
    for pattern in phospho_mods[3:6]:
        # Already counted if Phospho notation present
        if "(Phospho)" not in sequence:
            counts["pSer"] += len(re.findall(phospho_mods[3], sequence))
            break

    # Simpler approach: parse using known modification patterns
    counts = {"pSer": 0, "pThr": 0, "pTyr": 0}
    ser_patterns = [r"S\(Phospho\)", r"S\[80\]", r"S\[79\.966\]", r"S\[167\]"]
    thr_patterns = [r"T\(Phospho\)", r"T\[80\]", r"T\[79\.966\]", r"T\[181\]"]
    tyr_patterns = [r"Y\(Phospho\)", r"Y\[80\]", r"Y\[79\.966\]", r"Y\[243\]"]

    for pat in ser_patterns:
        counts["pSer"] += len(re.findall(pat, sequence))
    for pat in thr_patterns:
        counts["pThr"] += len(re.findall(pat, sequence))
    for pat in tyr_patterns:
        counts["pTyr"] += len(re.findall(pat, sequence))

    return counts


def get_peptide_length(sequence: str) -> int:
    """Get the amino acid length of a peptide using AASequence.

    Parameters
    ----------
    sequence:
        Peptide sequence string.

    Returns
    -------
    int
        Number of amino acid residues.
    """
    try:
        aa = oms.AASequence.fromString(sequence)
        return aa.size()
    except Exception:
        # Strip modifications manually for length estimation
        clean = re.sub(r"\[.*?\]", "", sequence)
        clean = re.sub(r"\(.*?\)", "", clean)
        return len(clean)


def compute_enrichment_stats(
    rows: List[Dict[str, str]],
) -> Tuple[Dict[str, int], Dict[str, float]]:
    """Compute enrichment statistics from search result rows.

    Parameters
    ----------
    rows:
        List of dicts with at least a 'sequence' key containing the peptide sequence.

    Returns
    -------
    tuple
        (counts, ratios) where counts has total, phospho, pSer, pThr, pTyr
        and ratios has enrichment_efficiency, pSer_ratio, pThr_ratio, pTyr_ratio.
    """
    counts = {"total": 0, "phospho": 0, "pSer": 0, "pThr": 0, "pTyr": 0}

    for row in rows:
        seq = row.get("sequence", row.get("peptide", ""))
        counts["total"] += 1
        if is_phosphopeptide(seq):
            counts["phospho"] += 1
            residue_counts = count_phospho_residues(seq)
            counts["pSer"] += residue_counts["pSer"]
            counts["pThr"] += residue_counts["pThr"]
            counts["pTyr"] += residue_counts["pTyr"]

    total_phospho_residues = counts["pSer"] + counts["pThr"] + counts["pTyr"]
    ratios: Dict[str, float] = {
        "enrichment_efficiency": counts["phospho"] / counts["total"] if counts["total"] > 0 else 0.0,
        "pSer_ratio": counts["pSer"] / total_phospho_residues if total_phospho_residues > 0 else 0.0,
        "pThr_ratio": counts["pThr"] / total_phospho_residues if total_phospho_residues > 0 else 0.0,
        "pTyr_ratio": counts["pTyr"] / total_phospho_residues if total_phospho_residues > 0 else 0.0,
    }

    return counts, ratios


def read_input(input_path: str) -> List[Dict[str, str]]:
    """Read search results TSV file."""
    rows = []
    with open(input_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def write_output(output_path: str, counts: Dict[str, int], ratios: Dict[str, float]) -> None:
    """Write enrichment report to TSV file."""
    with open(output_path, "w", newline="") as f:
        f.write("metric\tvalue\n")
        f.write(f"total_peptides\t{counts['total']}\n")
        f.write(f"phospho_peptides\t{counts['phospho']}\n")
        f.write(f"pSer_sites\t{counts['pSer']}\n")
        f.write(f"pThr_sites\t{counts['pThr']}\n")
        f.write(f"pTyr_sites\t{counts['pTyr']}\n")
        f.write(f"enrichment_efficiency\t{ratios['enrichment_efficiency']:.4f}\n")
        f.write(f"pSer_ratio\t{ratios['pSer_ratio']:.4f}\n")
        f.write(f"pThr_ratio\t{ratios['pThr_ratio']:.4f}\n")
        f.write(f"pTyr_ratio\t{ratios['pTyr_ratio']:.4f}\n")


@click.command(help="Compute phospho-enrichment efficiency and pSer/pThr/pTyr ratios.")
@click.option("--input", "input", required=True, help="Input search results TSV file")
@click.option("--output", required=True, help="Output enrichment report TSV file")
def main(input, output):
    rows = read_input(input)
    counts, ratios = compute_enrichment_stats(rows)
    write_output(output, counts, ratios)

    print(f"Total peptides:       {counts['total']}")
    print(f"Phospho peptides:     {counts['phospho']}")
    print(f"Enrichment efficiency: {ratios['enrichment_efficiency']:.4f}")
    print(f"pSer: {counts['pSer']}  pThr: {counts['pThr']}  pTyr: {counts['pTyr']}")
    print(f"pSer ratio: {ratios['pSer_ratio']:.4f}")
    print(f"pThr ratio: {ratios['pThr_ratio']:.4f}")
    print(f"pTyr ratio: {ratios['pTyr_ratio']:.4f}")


if __name__ == "__main__":
    main()
