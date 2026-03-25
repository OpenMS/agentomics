"""
RNA Digest
==========
In silico RNA digestion with common RNases.

Supported enzymes:
- RNase_T1: cleaves after G
- RNase_A: cleaves after C and U (pyrimidines)
- RNase_T2: cleaves after any nucleotide
- Cusativin: cleaves after C

Usage
-----
    python rna_digest.py --sequence AAUGCAAUGG --enzyme RNase_T1
    python rna_digest.py --sequence AAUGCAAUGG --enzyme RNase_A --missed-cleavages 1 --output fragments.tsv
"""

import argparse
import csv
import sys

try:
    import pyopenms as oms  # noqa: F401
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

PROTON = 1.007276

# Cleavage rules: enzyme -> set of nucleotides after which to cleave
ENZYME_RULES = {
    "RNase_T1": {"G"},
    "RNase_A": {"C", "U"},
    "RNase_T2": {"A", "C", "G", "U"},
    "Cusativin": {"C"},
}

# Monoisotopic residue masses for manual mass calculation
NUCLEOTIDE_RESIDUE_MASSES = {
    "A": 329.05252,
    "C": 305.04188,
    "G": 345.04744,
    "U": 306.02530,
}
WATER_MASS = 18.01056


def digest_rna(sequence: str, enzyme: str, missed_cleavages: int = 0) -> list:
    """Digest an RNA sequence in silico using the specified RNase.

    Parameters
    ----------
    sequence:
        RNA sequence (A, C, G, U).
    enzyme:
        Enzyme name (RNase_T1, RNase_A, RNase_T2, Cusativin).
    missed_cleavages:
        Maximum number of missed cleavages allowed.

    Returns
    -------
    list
        List of dicts with keys: fragment, start, end, missed_cleavages, mass.
    """
    sequence = sequence.upper().strip()
    if enzyme not in ENZYME_RULES:
        raise ValueError(f"Unknown enzyme: '{enzyme}'. Supported: {list(ENZYME_RULES.keys())}")

    for ch in sequence:
        if ch not in NUCLEOTIDE_RESIDUE_MASSES:
            raise ValueError(f"Invalid RNA nucleotide: '{ch}'.")

    cleavage_sites = ENZYME_RULES[enzyme]

    # Find cleavage positions (after these positions)
    cut_positions = []
    for i, nt in enumerate(sequence):
        if nt in cleavage_sites and i < len(sequence) - 1:
            cut_positions.append(i + 1)

    # Build basic fragments (0 missed cleavages)
    boundaries = [0] + cut_positions + [len(sequence)]

    fragments = []
    for mc in range(missed_cleavages + 1):
        for i in range(len(boundaries) - 1 - mc):
            start = boundaries[i]
            end = boundaries[i + 1 + mc]
            frag_seq = sequence[start:end]
            mass = _calculate_fragment_mass(frag_seq)
            fragments.append({
                "fragment": frag_seq,
                "start": start + 1,  # 1-based
                "end": end,
                "missed_cleavages": mc,
                "mass": mass,
            })

    return fragments


def _calculate_fragment_mass(fragment: str) -> float:
    """Calculate monoisotopic mass of an RNA fragment."""
    return sum(NUCLEOTIDE_RESIDUE_MASSES[nt] for nt in fragment) + WATER_MASS


def main():
    parser = argparse.ArgumentParser(description="In silico RNA digestion with RNases.")
    parser.add_argument("--sequence", required=True, help="RNA sequence (e.g. AAUGCAAUGG)")
    parser.add_argument(
        "--enzyme", required=True,
        choices=list(ENZYME_RULES.keys()),
        help="RNase enzyme name"
    )
    parser.add_argument("--missed-cleavages", type=int, default=0, help="Max missed cleavages (default: 0)")
    parser.add_argument("--output", help="Output TSV file (optional)")
    args = parser.parse_args()

    fragments = digest_rna(args.sequence, args.enzyme, args.missed_cleavages)

    if args.output:
        with open(args.output, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=["fragment", "start", "end", "missed_cleavages", "mass"],
                                    delimiter="\t")
            writer.writeheader()
            writer.writerows(fragments)
        print(f"Wrote {len(fragments)} fragments to {args.output}")
    else:
        print(f"Enzyme: {args.enzyme}")
        print(f"Sequence: {args.sequence}")
        print(f"Missed cleavages: {args.missed_cleavages}")
        print(f"\n{'Fragment':<20} {'Start':>5} {'End':>5} {'MC':>3} {'Mass':>12}")
        print("-" * 50)
        for f in fragments:
            print(f"{f['fragment']:<20} {f['start']:>5} {f['end']:>5} {f['missed_cleavages']:>3} "
                  f"{f['mass']:>12.4f}")


if __name__ == "__main__":
    main()
