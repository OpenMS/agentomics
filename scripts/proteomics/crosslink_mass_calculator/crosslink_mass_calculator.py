"""
Crosslink Mass Calculator
=========================
Calculate masses for crosslinked peptide pairs using common crosslinkers
(DSS, BS3, DSSO).

Supports:
- Built-in crosslinker mass table (DSS=138.068, BS3=138.068, DSSO=158.004)
- Multiple charge states
- TSV output

Usage
-----
    python crosslink_mass_calculator.py --peptide1 PEPTIDEK --peptide2 ANOTHERPEPTIDER --crosslinker DSS --charge 3
    python crosslink_mass_calculator.py --peptide1 PEPTIDEK --peptide2 ANOTHERPEPTIDER \\
        --crosslinker DSSO --output masses.tsv
"""

import argparse
import csv
import sys
from typing import Optional

try:
    import pyopenms as oms
except ImportError:
    sys.exit(
        "pyopenms is required. Install it with:  pip install pyopenms"
    )

PROTON = 1.007276

# Built-in crosslinker mass table (Da)
CROSSLINKER_MASSES = {
    "DSS": 138.068,
    "BS3": 138.068,
    "DSSO": 158.004,
}


def crosslinked_mass(
    peptide1: str,
    peptide2: str,
    crosslinker: str,
    charge: int = 1,
    custom_mass: Optional[float] = None,
) -> dict:
    """Calculate mass of a crosslinked peptide pair.

    Parameters
    ----------
    peptide1:
        First peptide sequence.
    peptide2:
        Second peptide sequence.
    crosslinker:
        Crosslinker name (DSS, BS3, DSSO) or custom name if custom_mass provided.
    charge:
        Charge state for m/z calculation.
    custom_mass:
        Optional custom crosslinker mass in Da.

    Returns
    -------
    dict
        Dictionary with peptide masses, crosslinker mass, total mass, and m/z.
    """
    if custom_mass is not None:
        xl_mass = custom_mass
    elif crosslinker.upper() in CROSSLINKER_MASSES:
        xl_mass = CROSSLINKER_MASSES[crosslinker.upper()]
    else:
        raise ValueError(
            f"Unknown crosslinker '{crosslinker}'. "
            f"Known: {', '.join(CROSSLINKER_MASSES.keys())}. "
            f"Provide --custom-mass for custom crosslinkers."
        )

    seq1 = oms.AASequence.fromString(peptide1)
    seq2 = oms.AASequence.fromString(peptide2)

    mass1 = seq1.getMonoWeight()
    mass2 = seq2.getMonoWeight()

    # Crosslinking forms a bond releasing water (two NH groups react with crosslinker)
    # Total mass = mass1 + mass2 + crosslinker_mass - 2*H2O
    # Actually, the standard model: XL mass = pep1 + pep2 + linker - 2*H2O is for
    # some chemistries. For NHS-ester crosslinkers (DSS, BS3, DSSO),
    # the reaction releases no extra water beyond what is already in the linker mass.
    # The crosslinker mass listed is the bridge mass after reaction.
    total_mass = mass1 + mass2 + xl_mass

    mz = (total_mass + charge * PROTON) / charge

    return {
        "peptide1": peptide1,
        "peptide2": peptide2,
        "crosslinker": crosslinker.upper() if custom_mass is None else crosslinker,
        "mass_peptide1": mass1,
        "mass_peptide2": mass2,
        "crosslinker_mass": xl_mass,
        "total_mass": total_mass,
        "charge": charge,
        "mz": mz,
    }


def write_tsv(results: list, output_path: str) -> None:
    """Write results to a TSV file.

    Parameters
    ----------
    results:
        List of result dicts from crosslinked_mass().
    output_path:
        Output file path.
    """
    if not results:
        return
    fieldnames = list(results[0].keys())
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in results:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(
        description="Calculate masses for crosslinked peptide pairs."
    )
    parser.add_argument("--peptide1", required=True, help="First peptide sequence")
    parser.add_argument("--peptide2", required=True, help="Second peptide sequence")
    parser.add_argument(
        "--crosslinker", required=True,
        help="Crosslinker name (DSS, BS3, DSSO) or custom name with --custom-mass"
    )
    parser.add_argument("--charge", type=int, default=1, help="Charge state (default: 1)")
    parser.add_argument("--custom-mass", type=float, default=None, help="Custom crosslinker mass in Da")
    parser.add_argument("--output", default=None, help="Output TSV file path")
    args = parser.parse_args()

    result = crosslinked_mass(
        args.peptide1, args.peptide2, args.crosslinker,
        charge=args.charge, custom_mass=args.custom_mass,
    )

    print(f"Peptide 1         : {result['peptide1']} ({result['mass_peptide1']:.6f} Da)")
    print(f"Peptide 2         : {result['peptide2']} ({result['mass_peptide2']:.6f} Da)")
    print(f"Crosslinker       : {result['crosslinker']} ({result['crosslinker_mass']:.6f} Da)")
    print(f"Total mass        : {result['total_mass']:.6f} Da")
    print(f"Charge            : {result['charge']}+")
    print(f"m/z               : {result['mz']:.6f}")

    if args.output:
        write_tsv([result], args.output)
        print(f"\nResults written to {args.output}")


if __name__ == "__main__":
    main()
