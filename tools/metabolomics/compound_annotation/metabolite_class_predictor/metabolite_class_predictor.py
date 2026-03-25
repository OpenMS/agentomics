"""
Metabolite Class Predictor
==========================
Predict compound class from molecular formula using mass defect, element
ratios (H:C, O:C), and Ring and Double Bond Equivalents (RDBE).

Uses heuristic classification rules based on:
- Mass defect ranges
- H:C and O:C ratios (van Krevelen-style)
- RDBE values
- Nitrogen/sulfur/phosphorus content

Usage
-----
    python metabolite_class_predictor.py --input formulas.tsv --output predictions.tsv
"""

import csv
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def get_element_counts(formula: str) -> dict[str, int]:
    """Extract element counts from a molecular formula using pyopenms.

    Parameters
    ----------
    formula:
        Empirical formula string, e.g. ``"C6H12O6"``.

    Returns
    -------
    dict
        Mapping element symbol to count, e.g. {"C": 6, "H": 12, "O": 6}.
    """
    ef = oms.EmpiricalFormula(formula)
    element_db = oms.ElementDB()
    counts = {}
    for element_name in ["Carbon", "Hydrogen", "Oxygen", "Nitrogen", "Sulfur", "Phosphorus"]:
        element = element_db.getElement(element_name)
        symbol = element.getSymbol()
        count = ef.getNumberOf(element)
        if count > 0:
            counts[symbol] = count
    return counts


def compute_exact_mass(formula: str) -> float:
    """Get the monoisotopic mass from formula.

    Parameters
    ----------
    formula:
        Empirical formula string.

    Returns
    -------
    float
        Monoisotopic mass.
    """
    ef = oms.EmpiricalFormula(formula)
    return ef.getMonoWeight()


def compute_mass_defect(mass: float) -> float:
    """Calculate mass defect (fractional part of the mass).

    Parameters
    ----------
    mass:
        Monoisotopic mass.

    Returns
    -------
    float
        Mass defect (value between 0 and 1).
    """
    return mass - int(mass)


def compute_rdbe(counts: dict[str, int]) -> float:
    """Calculate Ring and Double Bond Equivalents (RDBE).

    RDBE = 1 + C - H/2 + N/2

    Parameters
    ----------
    counts:
        Element counts dict.

    Returns
    -------
    float
        RDBE value.
    """
    c = counts.get("C", 0)
    h = counts.get("H", 0)
    n = counts.get("N", 0)
    return 1.0 + c - h / 2.0 + n / 2.0


def compute_element_ratios(counts: dict[str, int]) -> dict[str, float | None]:
    """Calculate H:C and O:C ratios.

    Parameters
    ----------
    counts:
        Element counts dict.

    Returns
    -------
    dict
        Keys: hc_ratio, oc_ratio. None if C == 0.
    """
    c = counts.get("C", 0)
    h = counts.get("H", 0)
    o = counts.get("O", 0)
    if c == 0:
        return {"hc_ratio": None, "oc_ratio": None}
    return {
        "hc_ratio": h / c,
        "oc_ratio": o / c,
    }


def classify_metabolite(
    formula: str,
) -> dict:
    """Predict compound class from molecular formula using heuristic rules.

    Classification is based on van Krevelen-style element ratio analysis:
    - Lipids: H:C > 1.5, O:C < 0.3, RDBE low
    - Carbohydrates: O:C 0.6-1.2, H:C 1.5-2.5
    - Amino acids / peptides: contains N, H:C 1.0-2.0, O:C 0.3-0.8
    - Nucleotides: contains N and P, O:C > 0.5
    - Terpenoids: H:C 1.0-1.8, O:C < 0.3, no N
    - Phenolics / polyketides: H:C 0.5-1.5, O:C 0.2-0.8, higher RDBE
    - Alkaloids: contains N, RDBE > 4, H:C < 1.5

    Parameters
    ----------
    formula:
        Empirical formula string.

    Returns
    -------
    dict
        Prediction results with keys: formula, exact_mass, mass_defect,
        rdbe, hc_ratio, oc_ratio, predicted_class, confidence.
    """
    counts = get_element_counts(formula)
    mass = compute_exact_mass(formula)
    mass_defect = compute_mass_defect(mass)
    rdbe = compute_rdbe(counts)
    ratios = compute_element_ratios(counts)
    hc = ratios["hc_ratio"]
    oc = ratios["oc_ratio"]

    has_n = counts.get("N", 0) > 0
    has_s = counts.get("S", 0) > 0
    has_p = counts.get("P", 0) > 0

    predicted_class = "Unknown"
    confidence = "low"

    if hc is not None and oc is not None:
        # Nucleotides: N + P, high O:C
        if has_n and has_p and oc > 0.5:
            predicted_class = "Nucleotide"
            confidence = "medium"

        # Carbohydrates: high O:C, high H:C, no N
        elif not has_n and 0.6 <= oc <= 1.2 and 1.5 <= hc <= 2.5:
            predicted_class = "Carbohydrate"
            confidence = "high" if 0.8 <= oc <= 1.1 and 1.8 <= hc <= 2.2 else "medium"

        # Amino acids / peptides: N present, moderate ratios
        elif has_n and 1.0 <= hc <= 2.2 and 0.2 <= oc <= 0.8 and rdbe < 6:
            predicted_class = "Amino acid / Peptide"
            confidence = "medium"

        # Alkaloids: N present, high aromaticity
        elif has_n and rdbe > 4 and hc < 1.5:
            predicted_class = "Alkaloid"
            confidence = "medium"

        # Lipids: high H:C, low O:C
        elif hc > 1.5 and oc < 0.3 and not has_n:
            predicted_class = "Lipid"
            confidence = "high" if hc > 1.7 and oc < 0.2 else "medium"

        # Terpenoids: moderate H:C, low O:C, no N
        elif 1.0 <= hc <= 1.8 and oc < 0.3 and not has_n:
            predicted_class = "Terpenoid"
            confidence = "medium"

        # Phenolics / polyketides: moderate ratios, higher RDBE
        elif 0.5 <= hc <= 1.5 and 0.2 <= oc <= 0.8 and rdbe > 3:
            predicted_class = "Phenolic / Polyketide"
            confidence = "medium"

        # Sulfur-containing organics
        elif has_s and has_n:
            predicted_class = "Sulfur-containing organic"
            confidence = "low"

    return {
        "formula": formula,
        "exact_mass": round(mass, 6),
        "mass_defect": round(mass_defect, 6),
        "rdbe": round(rdbe, 1),
        "hc_ratio": round(hc, 4) if hc is not None else None,
        "oc_ratio": round(oc, 4) if oc is not None else None,
        "has_nitrogen": has_n,
        "has_sulfur": has_s,
        "has_phosphorus": has_p,
        "predicted_class": predicted_class,
        "confidence": confidence,
    }


def classify_batch(formulas: list[dict]) -> list[dict]:
    """Classify a batch of formulas.

    Parameters
    ----------
    formulas:
        List of dicts with a 'formula' key.

    Returns
    -------
    list of dict
        Classification results.
    """
    results = []
    for row in formulas:
        formula = row.get("formula", "").strip()
        if not formula:
            continue
        result = classify_metabolite(formula)
        results.append(result)
    return results


def load_formulas(path: str) -> list[dict]:
    """Load formulas from TSV file. Expects a 'formula' column."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def write_predictions(predictions: list[dict], path: str) -> None:
    """Write predictions to TSV."""
    if not predictions:
        with open(path, "w") as fh:
            fh.write("# No predictions\n")
        return
    fieldnames = list(predictions[0].keys())
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(predictions)


@click.command()
@click.option("--input", "input_file", required=True, help="Formulas table (TSV) with 'formula' column")
@click.option("--output", required=True, help="Output predictions (TSV)")
def main(input_file, output) -> None:
    """CLI entry point."""
    formulas = load_formulas(input_file)
    predictions = classify_batch(formulas)
    write_predictions(predictions, output)
    print(f"Classified {len(predictions)} formulas, written to {output}")


if __name__ == "__main__":
    main()
