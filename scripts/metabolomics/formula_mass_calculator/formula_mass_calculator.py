"""
Formula Mass Calculator
========================
Calculate monoisotopic and average masses for molecular formulas,
with optional adduct correction. Supports batch mode via TSV input.

Usage
-----
    python formula_mass_calculator.py --formula C6H12O6 --adduct "[M+H]+" --output mass.json
    python formula_mass_calculator.py --batch formulas.tsv --output masses.tsv
"""

import argparse
import csv
import json
import sys

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

PROTON = 1.007276

# Common adduct definitions: name -> (mass_change, charge)
ADDUCT_TABLE = {
    "[M+H]+": (PROTON, 1),
    "[M+Na]+": (22.989218, 1),
    "[M+K]+": (38.963158, 1),
    "[M+NH4]+": (18.034164, 1),
    "[M-H]-": (-PROTON, 1),
    "[M+Cl]-": (34.969402, 1),
    "[M+FA-H]-": (44.998201, 1),
    "[M+2H]2+": (PROTON, 2),
    "[M-2H]2-": (-PROTON, 2),
    "[M]": (0.0, 1),
}


def calculate_formula_mass(formula: str, adduct: str = "[M]") -> dict:
    """Calculate masses for a molecular formula with an optional adduct.

    Parameters
    ----------
    formula:
        Molecular formula string, e.g. ``"C6H12O6"``.
    adduct:
        Adduct type, e.g. ``"[M+H]+"``.

    Returns
    -------
    dict
        Contains: formula, adduct, monoisotopic_mass, average_mass,
        mz (adduct-corrected m/z).
    """
    ef = oms.EmpiricalFormula(formula)
    mono = ef.getMonoWeight()
    avg = ef.getAverageWeight()

    if adduct in ADDUCT_TABLE:
        mass_add, charge = ADDUCT_TABLE[adduct]
    else:
        mass_add, charge = 0.0, 1

    if charge == 0:
        charge = 1

    mz = (mono + mass_add) / charge

    return {
        "formula": formula,
        "adduct": adduct,
        "monoisotopic_mass": round(mono, 6),
        "average_mass": round(avg, 6),
        "mz": round(mz, 6),
        "charge": charge,
    }


def batch_calculate(rows: list[dict]) -> list[dict]:
    """Calculate masses for a batch of formulas.

    Parameters
    ----------
    rows:
        List of dicts with keys: formula, and optionally adduct.

    Returns
    -------
    list[dict]
        Each dict is the output of ``calculate_formula_mass``.
    """
    results = []
    for row in rows:
        formula = row["formula"]
        adduct = row.get("adduct", "[M]").strip()
        if not adduct:
            adduct = "[M]"
        results.append(calculate_formula_mass(formula, adduct))
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Calculate masses for molecular formulas with adducts."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--formula", help="Single molecular formula (e.g. C6H12O6)")
    group.add_argument("--batch", metavar="FILE", help="Batch TSV file with formula column")
    parser.add_argument("--adduct", default="[M]", help='Adduct type (default: "[M]" = neutral)')
    parser.add_argument("--output", required=True, metavar="FILE", help="Output JSON or TSV file")
    args = parser.parse_args()

    if args.formula:
        result = calculate_formula_mass(args.formula, args.adduct)
        with open(args.output, "w") as fh:
            json.dump(result, fh, indent=2)
        print(f"Formula: {result['formula']}")
        print(f"Adduct:  {result['adduct']}")
        print(f"Mono:    {result['monoisotopic_mass']:.6f} Da")
        print(f"Avg:     {result['average_mass']:.6f} Da")
        print(f"m/z:     {result['mz']:.6f}")
    else:
        rows = []
        with open(args.batch) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                rows.append(row)

        results = batch_calculate(rows)

        with open(args.output, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=["formula", "adduct", "monoisotopic_mass", "average_mass", "mz", "charge"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(results)

        print(f"Calculated masses for {len(results)} formulas, written to {args.output}")


if __name__ == "__main__":
    main()
