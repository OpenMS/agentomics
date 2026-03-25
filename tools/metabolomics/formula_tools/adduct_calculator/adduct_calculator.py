"""
Adduct Calculator
=================
Compute m/z values for all common ESI adducts given a molecular formula
or neutral monoisotopic mass.

Includes built-in tables for positive-mode and negative-mode adducts.

Usage
-----
    python adduct_calculator.py --formula C6H12O6 --mode positive --output adducts.tsv
    python adduct_calculator.py --mass 180.0634 --mode negative --output adducts.tsv
"""

import csv
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

PROTON = 1.007276
ELECTRON = 0.000549

# Adduct definitions: (name, mass_change, charge_multiplier)
POSITIVE_ADDUCTS = [
    ("[M+H]+", PROTON, 1),
    ("[M+Na]+", 22.989218, 1),
    ("[M+K]+", 38.963158, 1),
    ("[M+NH4]+", 18.034164, 1),
    ("[M+2H]2+", PROTON, 2),
    ("[M+H+Na]2+", (PROTON + 22.989218) / 1.0, 2),
    ("[M+Li]+", 7.016003, 1),
    ("[M+CH3OH+H]+", 33.033491, 1),
    ("[M+ACN+H]+", 42.033823, 1),
    ("[M+2Na-H]+", 2 * 22.989218 - PROTON, 1),
]

NEGATIVE_ADDUCTS = [
    ("[M-H]-", -PROTON, 1),
    ("[M+Cl]-", 34.969402, 1),
    ("[M+FA-H]-", 44.998201, 1),
    ("[M+CH3COO]-", 59.013851, 1),
    ("[M-2H]2-", -PROTON, 2),
    ("[M+Br]-", 78.918885, 1),
    ("[M+Na-2H]-", 22.989218 - 2 * PROTON, 1),
]


def compute_adducts(neutral_mass: float, mode: str = "positive") -> list[dict]:
    """Compute m/z for all adducts of the given neutral mass.

    Parameters
    ----------
    neutral_mass:
        Monoisotopic neutral mass in Da.
    mode:
        ``"positive"`` or ``"negative"``.

    Returns
    -------
    list[dict]
        Each dict has: adduct, mz, charge.
    """
    adducts = POSITIVE_ADDUCTS if mode == "positive" else NEGATIVE_ADDUCTS
    results = []

    for name, mass_add, charge in adducts:
        mz = (neutral_mass + mass_add) / charge
        results.append({
            "adduct": name,
            "mz": round(mz, 6),
            "charge": charge,
        })

    return results


def formula_to_mass(formula: str) -> float:
    """Convert a molecular formula string to monoisotopic mass using pyopenms.

    Parameters
    ----------
    formula:
        Molecular formula, e.g. ``"C6H12O6"``.

    Returns
    -------
    float
        Monoisotopic mass in Da.
    """
    ef = oms.EmpiricalFormula(formula)
    return ef.getMonoWeight()


@click.command()
@click.option("--formula", default=None, help="Molecular formula (e.g. C6H12O6)")
@click.option("--mass", type=float, default=None, help="Neutral monoisotopic mass in Da")
@click.option("--mode", type=click.Choice(["positive", "negative"]), default="positive",
              help="Ionization mode (default: positive)")
@click.option("--output", required=True, help="Output TSV file")
def main(formula, mass, mode, output):
    if not formula and mass is None:
        raise click.UsageError("Either --formula or --mass must be provided.")
    if formula and mass is not None:
        raise click.UsageError("--formula and --mass are mutually exclusive.")

    if formula:
        mass = formula_to_mass(formula)
        print(f"Formula: {formula}  Mass: {mass:.6f} Da")
    else:
        print(f"Mass: {mass:.6f} Da")

    adducts = compute_adducts(mass, mode)

    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["adduct", "mz", "charge"], delimiter="\t")
        writer.writeheader()
        writer.writerows(adducts)

    print(f"\n{len(adducts)} adducts written to {output}")
    for a in adducts:
        print(f"  {a['adduct']:<20}  m/z = {a['mz']:.6f}  (z={a['charge']})")


if __name__ == "__main__":
    main()
