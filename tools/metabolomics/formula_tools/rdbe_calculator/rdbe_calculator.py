"""
RDBE Calculator
================
Calculate Ring and Double Bond Equivalents (RDBE) for molecular formulas.
RDBE = (2C + 2 - H + N + P) / 2

Usage
-----
    python rdbe_calculator.py --input formulas.tsv --output rdbe.tsv
"""

import csv
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def get_element_counts(formula: str) -> dict:
    """Extract element counts from a molecular formula.

    Parameters
    ----------
    formula:
        Molecular formula string.

    Returns
    -------
    dict mapping element symbols (str) to counts (int).
    """
    ef = oms.EmpiricalFormula(formula)
    composition = ef.getElementalComposition()
    return {k.decode(): v for k, v in composition.items()}


def calculate_rdbe(formula: str) -> float:
    """Calculate RDBE for a molecular formula.

    RDBE = (2C + 2 - H + N + P) / 2

    Parameters
    ----------
    formula:
        Molecular formula string, e.g. ``"C6H6"``.

    Returns
    -------
    float: RDBE value.
    """
    counts = get_element_counts(formula)
    c = counts.get("C", 0)
    h = counts.get("H", 0)
    n = counts.get("N", 0)
    p = counts.get("P", 0)
    return (2 * c + 2 - h + n + p) / 2.0


def calculate_rdbe_batch(formulas: list) -> list:
    """Calculate RDBE for a list of formulas.

    Parameters
    ----------
    formulas:
        List of molecular formula strings.

    Returns
    -------
    list of dicts with keys: formula, rdbe, C, H, N, P.
    """
    results = []
    for formula in formulas:
        counts = get_element_counts(formula)
        rdbe = calculate_rdbe(formula)
        results.append({
            "formula": formula,
            "C": counts.get("C", 0),
            "H": counts.get("H", 0),
            "N": counts.get("N", 0),
            "P": counts.get("P", 0),
            "rdbe": round(rdbe, 1),
        })
    return results


@click.command()
@click.option("--input", "input_file", required=True, help="TSV file with a 'formula' column.")
@click.option("--output", required=True, help="Output TSV file with RDBE values.")
def main(input_file, output) -> None:
    """CLI entry point."""
    formulas = []
    with open(input_file, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            formulas.append(row["formula"])

    results = calculate_rdbe_batch(formulas)

    fieldnames = ["formula", "C", "H", "N", "P", "rdbe"]
    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    print(f"Calculated RDBE for {len(results)} formulas, wrote to {output}")


if __name__ == "__main__":
    main()
