"""
Immunopeptidome QC
==================
Quality control for immunopeptidomics data: peptide length distribution,
anchor residue frequencies, and per-position information content.

HLA class I peptides are expected to be 8-12 amino acids; HLA class II
peptides 12-25 amino acids.  The tool reads a TSV with a ``sequence`` column,
validates lengths against the expected range, computes positional amino acid
frequencies and information content (bits), and reports anchor residue
enrichment.

Usage
-----
    python immunopeptidome_qc.py --input hla_peptides.tsv --hla-class I \
        --output length_dist.tsv --motifs anchor_freq.tsv
"""

import csv
import math
import sys
from collections import Counter
from typing import Dict, List, Tuple

import click
import pyopenms as oms

# Expected peptide length ranges per HLA class
HLA_LENGTH_RANGES: Dict[str, Tuple[int, int]] = {
    "I": (8, 12),
    "II": (12, 25),
}

# Number of standard amino acids
NUM_AA = 20


def validate_sequence(sequence: str) -> bool:
    """Return True if *sequence* can be parsed by pyopenms AASequence."""
    try:
        aa = oms.AASequence.fromString(sequence)
        return aa.size() > 0
    except Exception:
        return False


def length_distribution(sequences: List[str]) -> Dict[int, int]:
    """Compute a histogram of peptide lengths.

    Parameters
    ----------
    sequences:
        List of amino acid sequences.

    Returns
    -------
    dict
        Mapping of peptide length to count, sorted by length.
    """
    counter: Counter = Counter()
    for seq in sequences:
        aa = oms.AASequence.fromString(seq)
        counter[aa.size()] += 1
    return dict(sorted(counter.items()))


def length_qc(length_dist: Dict[int, int], hla_class: str) -> Dict[str, object]:
    """Evaluate length distribution against expected HLA range.

    Returns
    -------
    dict
        ``in_range_count``, ``out_of_range_count``, ``in_range_fraction``,
        ``expected_min``, ``expected_max``.
    """
    lo, hi = HLA_LENGTH_RANGES[hla_class]
    in_range = sum(cnt for length, cnt in length_dist.items() if lo <= length <= hi)
    total = sum(length_dist.values())
    return {
        "expected_min": lo,
        "expected_max": hi,
        "in_range_count": in_range,
        "out_of_range_count": total - in_range,
        "in_range_fraction": in_range / total if total > 0 else 0.0,
    }


def positional_frequencies(sequences: List[str], target_length: int) -> List[Dict[str, float]]:
    """Compute per-position amino acid frequencies for peptides of *target_length*.

    Parameters
    ----------
    sequences:
        Input peptide sequences.
    target_length:
        Only peptides of exactly this length are considered.

    Returns
    -------
    list of dict
        One dict per position mapping residue one-letter code to frequency.
    """
    counters: List[Counter] = [Counter() for _ in range(target_length)]
    n = 0
    for seq in sequences:
        aa = oms.AASequence.fromString(seq)
        if aa.size() != target_length:
            continue
        n += 1
        for i in range(target_length):
            residue = aa.getResidue(i).getOneLetterCode()
            counters[i][residue] += 1
    if n == 0:
        return [{} for _ in range(target_length)]
    result: List[Dict[str, float]] = []
    for counter in counters:
        total = sum(counter.values())
        result.append({aa: cnt / total for aa, cnt in sorted(counter.items())})
    return result


def information_content(freq: Dict[str, float]) -> float:
    """Shannon information content in bits for a single position.

    IC = log2(N) + sum(p * log2(p))  where N = 20 amino acids.
    """
    max_entropy = math.log2(NUM_AA)
    entropy = 0.0
    for p in freq.values():
        if p > 0:
            entropy -= p * math.log2(p)
    return max_entropy - entropy


def anchor_residue_frequencies(
    sequences: List[str], hla_class: str
) -> Dict[str, Dict[str, float]]:
    """Compute anchor residue frequencies for dominant peptide lengths.

    For HLA-I the canonical anchors are position 2 and the C-terminal position
    (for 9-mers).  For HLA-II the anchors are positions 1, 4, 6, 9 of the
    9-mer binding core, but since core identification is non-trivial we report
    position 1 and the C-terminal position for the most common length.

    Returns
    -------
    dict
        Keys are position labels (e.g. ``"P2"``, ``"PC"``), values are dicts
        mapping residue to frequency.
    """
    # find most common length within range
    lo, hi = HLA_LENGTH_RANGES[hla_class]
    length_counter: Counter = Counter()
    for seq in sequences:
        aa = oms.AASequence.fromString(seq)
        sz = aa.size()
        if lo <= sz <= hi:
            length_counter[sz] += 1
    if not length_counter:
        return {}
    dominant_length = length_counter.most_common(1)[0][0]

    freqs = positional_frequencies(sequences, dominant_length)
    anchors: Dict[str, Dict[str, float]] = {}
    if hla_class == "I":
        if len(freqs) >= 2:
            anchors["P2"] = freqs[1]
        anchors["PC"] = freqs[-1]
    else:
        anchors["P1"] = freqs[0]
        anchors["PC"] = freqs[-1]
    return anchors


def run_qc(
    sequences: List[str], hla_class: str
) -> Tuple[Dict[int, int], Dict[str, object], Dict[str, Dict[str, float]], List[float]]:
    """Run the full QC pipeline.

    Returns
    -------
    tuple
        (length_dist, length_qc_result, anchor_freqs, info_content_per_position)
    """
    dist = length_distribution(sequences)
    qc = length_qc(dist, hla_class)

    anchors = anchor_residue_frequencies(sequences, hla_class)

    # info content for dominant length
    lo, hi = HLA_LENGTH_RANGES[hla_class]
    length_counter: Counter = Counter()
    for seq in sequences:
        aa = oms.AASequence.fromString(seq)
        sz = aa.size()
        if lo <= sz <= hi:
            length_counter[sz] += 1
    if length_counter:
        dominant = length_counter.most_common(1)[0][0]
        freqs = positional_frequencies(sequences, dominant)
        ic_values = [information_content(f) for f in freqs]
    else:
        ic_values = []

    return dist, qc, anchors, ic_values


@click.command(help="QC for immunopeptidomics: length distribution, anchor residue frequencies, information content.")
@click.option("--input", "input", required=True, help="Input TSV with 'sequence' column")
@click.option("--hla-class", required=True, type=click.Choice(["I", "II"]), help="HLA class (I or II)")
@click.option("--output", required=True, help="Output TSV for length distribution")
@click.option("--motifs", required=True, help="Output TSV for anchor residue frequencies")
def main(input, hla_class, output, motifs) -> None:
    # Read sequences
    sequences: List[str] = []
    with open(input, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            seq = row.get("sequence", "").strip()
            if seq and validate_sequence(seq):
                sequences.append(seq)

    if not sequences:
        sys.exit("No valid sequences found in input file.")

    dist, qc, anchors, ic_values = run_qc(sequences, hla_class)

    # Write length distribution
    with open(output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["length", "count"])
        for length, count in dist.items():
            writer.writerow([length, count])
        writer.writerow([])
        writer.writerow(["metric", "value"])
        for key, val in qc.items():
            writer.writerow([key, val])
        if ic_values:
            writer.writerow([])
            writer.writerow(["position", "information_content_bits"])
            for i, ic in enumerate(ic_values, 1):
                writer.writerow([i, f"{ic:.4f}"])

    # Write anchor frequencies
    with open(motifs, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["anchor_position", "residue", "frequency"])
        for pos_label, freq_dict in anchors.items():
            for residue, freq in sorted(freq_dict.items()):
                writer.writerow([pos_label, residue, f"{freq:.4f}"])

    print(f"Length distribution written to {output}")
    print(f"Anchor frequencies written to {motifs}")
    print(f"Total peptides: {sum(dist.values())}, in-range: {qc['in_range_fraction']:.1%}")


if __name__ == "__main__":
    main()
