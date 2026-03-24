"""
Mass Accuracy Calculator
========================
Calculate the mass accuracy (ppm error) between a theoretical mass
(or peptide/formula string) and one or more observed m/z values.

Supports both:
- Amino acid sequence inputs  (e.g. ``PEPTIDEK``)
- Empirical formula inputs    (e.g. ``C6H12O6``)

Usage
-----
    # Peptide sequence
    python mass_accuracy_calculator.py --sequence PEPTIDEK --observed 803.4560

    # Molecular formula (charge 1 default)
    python mass_accuracy_calculator.py --formula C6H12O6 --observed 181.0709

    # Multiple observed values at charge 2
    python mass_accuracy_calculator.py --sequence ACDEFGHIK --charge 2 \\
        --observed 554.2478 554.2480 554.2482
"""

import argparse
import sys

try:
    import pyopenms as oms
except ImportError:
    sys.exit(
        "pyopenms is required. Install it with:  pip install pyopenms"
    )

PROTON = 1.007276


def theoretical_mz_from_sequence(sequence: str, charge: int) -> float:
    """Compute the theoretical m/z for a peptide sequence.

    Parameters
    ----------
    sequence:
        Amino acid sequence, optionally with bracket-enclosed modifications.
    charge:
        Charge state.

    Returns
    -------
    float
        Theoretical m/z value.
    """
    aa_seq = oms.AASequence.fromString(sequence)
    mass = aa_seq.getMonoWeight()
    return (mass + charge * PROTON) / charge


def theoretical_mz_from_formula(formula: str, charge: int) -> float:
    """Compute the theoretical m/z for a molecular formula.

    Parameters
    ----------
    formula:
        Empirical formula string, e.g. ``"C6H12O6"``.
    charge:
        Charge state (used for proton addition).

    Returns
    -------
    float
        Theoretical m/z value (monoisotopic).
    """
    ef = oms.EmpiricalFormula(formula)
    mass = ef.getMonoWeight()
    return (mass + charge * PROTON) / charge


def ppm_error(theoretical: float, observed: float) -> float:
    """Calculate the mass accuracy in parts-per-million (ppm).

    Parameters
    ----------
    theoretical:
        Theoretical m/z value.
    observed:
        Observed m/z value.

    Returns
    -------
    float
        PPM error (positive = observed > theoretical).
    """
    return (observed - theoretical) / theoretical * 1e6


def main():
    parser = argparse.ArgumentParser(
        description="Compute m/z mass accuracy (ppm error) using pyopenms."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--sequence",
        help="Peptide sequence (e.g. PEPTIDEK)",
    )
    group.add_argument(
        "--formula",
        help="Molecular formula (e.g. C6H12O6)",
    )
    parser.add_argument(
        "--charge",
        type=int,
        default=1,
        help="Charge state (default: 1)",
    )
    parser.add_argument(
        "--observed",
        nargs="+",
        type=float,
        required=True,
        metavar="MZ",
        help="Observed m/z value(s)",
    )
    args = parser.parse_args()

    if args.sequence:
        theoretical = theoretical_mz_from_sequence(args.sequence, args.charge)
        label = f"sequence={args.sequence}"
    else:
        theoretical = theoretical_mz_from_formula(args.formula, args.charge)
        label = f"formula={args.formula}"

    print(f"Theoretical m/z ({label}, charge {args.charge}+): {theoretical:.6f}")
    print(f"\n{'Observed m/z':>14}  {'PPM error':>10}")
    print("-" * 28)
    for obs in args.observed:
        ppm = ppm_error(theoretical, obs)
        print(f"{obs:>14.6f}  {ppm:>+10.4f}")


if __name__ == "__main__":
    main()
