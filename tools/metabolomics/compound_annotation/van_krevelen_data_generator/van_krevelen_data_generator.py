"""
Van Krevelen Data Generator
============================
Compute H:C and O:C ratios from molecular formulas and classify compounds
into biochemical classes based on their position in Van Krevelen space.

Usage
-----
    python van_krevelen_data_generator.py --input formulas.tsv --classify --output van_krevelen.tsv
"""

import csv

import click
import pyopenms as oms

# Biochemical class regions defined by H:C and O:C ratio boundaries
BIOCHEMICAL_CLASSES = {
    "lipids": {"hc_min": 1.5, "hc_max": 2.5, "oc_min": 0.0, "oc_max": 0.3},
    "carbohydrates": {"hc_min": 1.5, "hc_max": 2.5, "oc_min": 0.6, "oc_max": 1.2},
    "amino_acids": {"hc_min": 1.4, "hc_max": 2.0, "oc_min": 0.3, "oc_max": 0.7},
    "nucleotides": {"hc_min": 1.0, "hc_max": 1.5, "oc_min": 0.5, "oc_max": 1.0},
}


def compute_ratios(formula: str) -> dict:
    """Compute H:C and O:C ratios from a molecular formula.

    Parameters
    ----------
    formula:
        Empirical formula string, e.g. ``"C6H12O6"``.

    Returns
    -------
    dict with keys: formula, C, H, O, hc_ratio, oc_ratio
    """
    ef = oms.EmpiricalFormula(formula)
    composition = ef.getElementalComposition()

    c_count = composition.get(b"C", 0)
    h_count = composition.get(b"H", 0)
    o_count = composition.get(b"O", 0)

    if c_count == 0:
        raise ValueError(f"Formula '{formula}' contains no carbon atoms; cannot compute ratios.")

    hc_ratio = h_count / c_count
    oc_ratio = o_count / c_count

    return {
        "formula": formula,
        "C": c_count,
        "H": h_count,
        "O": o_count,
        "hc_ratio": round(hc_ratio, 4),
        "oc_ratio": round(oc_ratio, 4),
    }


def classify_compound(hc_ratio: float, oc_ratio: float) -> str:
    """Classify a compound into a biochemical class based on Van Krevelen regions.

    Parameters
    ----------
    hc_ratio:
        Hydrogen-to-carbon ratio.
    oc_ratio:
        Oxygen-to-carbon ratio.

    Returns
    -------
    str: The biochemical class name, or "unclassified" if no region matches.
    """
    for cls_name, bounds in BIOCHEMICAL_CLASSES.items():
        if (bounds["hc_min"] <= hc_ratio <= bounds["hc_max"]
                and bounds["oc_min"] <= oc_ratio <= bounds["oc_max"]):
            return cls_name
    return "unclassified"


def process_formulas(formulas: list, classify: bool = False) -> list:
    """Process a list of formulas and optionally classify them.

    Parameters
    ----------
    formulas:
        List of molecular formula strings.
    classify:
        Whether to add biochemical class assignment.

    Returns
    -------
    list of dicts with ratio data and optional classification.
    """
    results = []
    for formula in formulas:
        row = compute_ratios(formula)
        if classify:
            row["class"] = classify_compound(row["hc_ratio"], row["oc_ratio"])
        results.append(row)
    return results


@click.command()
@click.option("--input", "input_file", required=True, help="TSV file with a 'formula' column.")
@click.option("--classify", is_flag=True, help="Add biochemical class assignment.")
@click.option("--output", required=True, help="Output TSV file with ratios.")
def main(input_file, classify, output) -> None:
    """CLI entry point."""
    formulas = []
    with open(input_file, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            formulas.append(row["formula"])

    results = process_formulas(formulas, classify=classify)

    fieldnames = ["formula", "C", "H", "O", "hc_ratio", "oc_ratio"]
    if classify:
        fieldnames.append("class")

    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    print(f"Wrote {len(results)} entries to {output}")


if __name__ == "__main__":
    main()
