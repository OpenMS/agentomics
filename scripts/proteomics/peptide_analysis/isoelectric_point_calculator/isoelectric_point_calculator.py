"""
Isoelectric Point Calculator
==============================
Calculate pI for peptides and proteins using Henderson-Hasselbalch equation.

Features
--------
- Multiple pKa sets: Lehninger, EMBOSS, Stryer, Solomon
- Bisection method for accurate pI determination
- Batch processing from FASTA or TSV input
- Charge curve generation at multiple pH values

Usage
-----
    python isoelectric_point_calculator.py --sequence ACDEFGHIK --pk-set lehninger --output pi.json
    python isoelectric_point_calculator.py --fasta proteins.fasta --output pi.tsv
"""

import argparse
import csv
import json
import sys

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

# pKa sets from different sources
PKA_SETS = {
    "lehninger": {
        "nterm": 9.69, "cterm": 2.34,
        "C": 8.33, "D": 3.65, "E": 4.25, "H": 6.00, "K": 10.53, "R": 12.48, "Y": 10.07,
    },
    "emboss": {
        "nterm": 8.6, "cterm": 3.6,
        "C": 8.5, "D": 3.9, "E": 4.1, "H": 6.5, "K": 10.8, "R": 12.5, "Y": 10.1,
    },
    "stryer": {
        "nterm": 8.0, "cterm": 3.1,
        "C": 8.3, "D": 4.0, "E": 4.1, "H": 6.0, "K": 10.0, "R": 12.5, "Y": 10.1,
    },
    "solomon": {
        "nterm": 9.56, "cterm": 2.34,
        "C": 8.00, "D": 3.65, "E": 4.25, "H": 6.00, "K": 10.53, "R": 12.48, "Y": 10.07,
    },
}

# Residues that contribute positive charge when protonated
POSITIVE_RESIDUES = {"H", "K", "R"}
# Residues that contribute negative charge when deprotonated
NEGATIVE_RESIDUES = {"C", "D", "E", "Y"}


def charge_at_ph(sequence: str, ph: float, pk_set: str = "lehninger") -> float:
    """Calculate the net charge of a peptide at a given pH.

    Parameters
    ----------
    sequence : str
        One-letter amino acid sequence.
    ph : float
        pH value.
    pk_set : str
        pKa set to use.

    Returns
    -------
    float
        Net charge at the given pH.
    """
    pka = PKA_SETS.get(pk_set, PKA_SETS["lehninger"])
    charge = 0.0

    # N-terminus (positive when protonated)
    charge += 1.0 / (1.0 + 10 ** (ph - pka["nterm"]))
    # C-terminus (negative when deprotonated)
    charge -= 1.0 / (1.0 + 10 ** (pka["cterm"] - ph))

    for aa in sequence:
        if aa in POSITIVE_RESIDUES:
            charge += 1.0 / (1.0 + 10 ** (ph - pka[aa]))
        elif aa in NEGATIVE_RESIDUES:
            charge -= 1.0 / (1.0 + 10 ** (pka[aa] - ph))

    return charge


def calculate_pi(sequence: str, pk_set: str = "lehninger", precision: float = 0.01) -> float:
    """Calculate the isoelectric point using bisection.

    Parameters
    ----------
    sequence : str
        Amino acid sequence (plain one-letter code).
    pk_set : str
        pKa set name.
    precision : float
        Desired precision.

    Returns
    -------
    float
        Estimated pI.
    """
    low, high = 0.0, 14.0
    while (high - low) > precision:
        mid = (low + high) / 2.0
        if charge_at_ph(sequence, mid, pk_set) > 0:
            low = mid
        else:
            high = mid
    return round((low + high) / 2.0, 2)


def calculate_charge_curve(sequence: str, pk_set: str = "lehninger",
                           ph_start: float = 0.0, ph_end: float = 14.0,
                           ph_step: float = 0.5) -> list:
    """Calculate charge across a range of pH values.

    Parameters
    ----------
    sequence : str
        Amino acid sequence.
    pk_set : str
        pKa set name.
    ph_start : float
        Starting pH.
    ph_end : float
        Ending pH.
    ph_step : float
        pH increment.

    Returns
    -------
    list
        List of (pH, charge) tuples.
    """
    curve = []
    ph = ph_start
    while ph <= ph_end + 0.001:
        c = charge_at_ph(sequence, ph, pk_set)
        curve.append({"ph": round(ph, 2), "charge": round(c, 4)})
        ph += ph_step
    return curve


def calculate_pi_from_sequence(sequence: str, pk_set: str = "lehninger") -> dict:
    """Calculate pI and related properties for a sequence.

    Parameters
    ----------
    sequence : str
        Peptide or protein sequence (plain or pyopenms notation).
    pk_set : str
        pKa set name.

    Returns
    -------
    dict
        Dictionary with pI, charge at pI, and sequence info.
    """
    aa_seq = oms.AASequence.fromString(sequence)
    plain = aa_seq.toUnmodifiedString()
    pi = calculate_pi(plain, pk_set)
    charge = charge_at_ph(plain, pi, pk_set)

    return {
        "sequence": sequence,
        "unmodified_sequence": plain,
        "length": len(plain),
        "pI": pi,
        "charge_at_pI": round(charge, 4),
        "pk_set": pk_set,
    }


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(description="Calculate isoelectric point for peptides/proteins.")
    parser.add_argument("--sequence", type=str, help="Single amino acid sequence.")
    parser.add_argument("--fasta", type=str, help="FASTA file with protein sequences.")
    parser.add_argument("--pk-set", choices=list(PKA_SETS.keys()), default="lehninger",
                        help="pKa value set (default: lehninger).")
    parser.add_argument("--charge-curve", action="store_true", help="Also output charge curve.")
    parser.add_argument("--output", type=str, help="Output file (.json or .tsv).")
    args = parser.parse_args()

    if not args.sequence and not args.fasta:
        parser.error("Provide --sequence or --fasta.")

    results = []
    if args.sequence:
        result = calculate_pi_from_sequence(args.sequence, args.pk_set)
        if args.charge_curve:
            aa_seq = oms.AASequence.fromString(args.sequence)
            result["charge_curve"] = calculate_charge_curve(aa_seq.toUnmodifiedString(), args.pk_set)
        results.append(result)
    elif args.fasta:
        entries = []
        oms.FASTAFile().load(args.fasta, entries)
        for entry in entries:
            result = calculate_pi_from_sequence(entry.sequence, args.pk_set)
            result["accession"] = entry.identifier
            results.append(result)

    if args.output:
        if args.output.endswith(".json"):
            with open(args.output, "w") as fh:
                json.dump(results if len(results) > 1 else results[0], fh, indent=2)
        else:
            with open(args.output, "w", newline="") as fh:
                fieldnames = ["sequence", "length", "pI", "charge_at_pI", "pk_set"]
                if args.fasta:
                    fieldnames.insert(0, "accession")
                writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
                writer.writeheader()
                writer.writerows(results)
        print(f"Results written to {args.output}")
    else:
        for r in results:
            acc = r.get("accession", r.get("sequence", ""))
            print(f"{acc}\tpI={r['pI']}\tcharge@pI={r['charge_at_pI']}")


if __name__ == "__main__":
    main()
