"""
Formula Validator - Seven Golden Rules
========================================
Apply the Seven Golden Rules (Kind & Fiehn, BMC Bioinformatics, 2007) to
filter molecular formulas for plausibility.

Usage
-----
    python formula_validator_golden_rules.py --input formulas.tsv --rules all --output validated.tsv
"""

import csv
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def get_element_counts(formula: str) -> dict:
    """Extract element counts from a molecular formula using pyopenms.

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


def compute_rdbe(counts: dict) -> float:
    """Compute the Ring and Double Bond Equivalents (RDBE).

    RDBE = (2C + 2 - H + N + P) / 2

    Parameters
    ----------
    counts:
        Element counts dict.

    Returns
    -------
    float: RDBE value.
    """
    c = counts.get("C", 0)
    h = counts.get("H", 0)
    n = counts.get("N", 0)
    p = counts.get("P", 0)
    return (2 * c + 2 - h + n + p) / 2.0


def check_rdbe_nonnegative(counts: dict) -> bool:
    """Rule: RDBE must be non-negative.

    Parameters
    ----------
    counts:
        Element counts dict.

    Returns
    -------
    bool: True if valid.
    """
    return compute_rdbe(counts) >= 0


def check_hc_ratio(counts: dict) -> bool:
    """Rule: H/C ratio must be between 0.2 and 3.1.

    Parameters
    ----------
    counts:
        Element counts dict.

    Returns
    -------
    bool: True if valid.
    """
    c = counts.get("C", 0)
    h = counts.get("H", 0)
    if c == 0:
        return False
    ratio = h / c
    return 0.2 <= ratio <= 3.1


def check_nc_ratio(counts: dict) -> bool:
    """Rule: N/C ratio must be <= 1.3.

    Parameters
    ----------
    counts:
        Element counts dict.

    Returns
    -------
    bool: True if valid.
    """
    c = counts.get("C", 0)
    n = counts.get("N", 0)
    if c == 0:
        return n == 0
    return n / c <= 1.3


def check_oc_ratio(counts: dict) -> bool:
    """Rule: O/C ratio must be <= 1.2.

    Parameters
    ----------
    counts:
        Element counts dict.

    Returns
    -------
    bool: True if valid.
    """
    c = counts.get("C", 0)
    o = counts.get("O", 0)
    if c == 0:
        return o == 0
    return o / c <= 1.2


def check_sc_ratio(counts: dict) -> bool:
    """Rule: S/C ratio must be <= 0.8.

    Parameters
    ----------
    counts:
        Element counts dict.

    Returns
    -------
    bool: True if valid.
    """
    c = counts.get("C", 0)
    s = counts.get("S", 0)
    if c == 0:
        return s == 0
    return s / c <= 0.8


def check_pc_ratio(counts: dict) -> bool:
    """Rule: P/C ratio must be <= 0.3.

    Parameters
    ----------
    counts:
        Element counts dict.

    Returns
    -------
    bool: True if valid.
    """
    c = counts.get("C", 0)
    p = counts.get("P", 0)
    if c == 0:
        return p == 0
    return p / c <= 0.3


# All available rules
ALL_RULES = {
    "rdbe": check_rdbe_nonnegative,
    "hc": check_hc_ratio,
    "nc": check_nc_ratio,
    "oc": check_oc_ratio,
    "sc": check_sc_ratio,
    "pc": check_pc_ratio,
}


def validate_formula(formula: str, rules: list) -> dict:
    """Validate a formula against the specified rules.

    Parameters
    ----------
    formula:
        Molecular formula string.
    rules:
        List of rule names to apply, or ``["all"]`` for all rules.

    Returns
    -------
    dict with formula, per-rule pass/fail, and overall valid flag.
    """
    counts = get_element_counts(formula)
    rdbe_val = compute_rdbe(counts)

    if "all" in rules:
        active_rules = list(ALL_RULES.keys())
    else:
        active_rules = [r for r in rules if r in ALL_RULES]

    result = {"formula": formula, "rdbe": round(rdbe_val, 1)}
    all_pass = True
    for rule_name in active_rules:
        passed = ALL_RULES[rule_name](counts)
        result[f"rule_{rule_name}"] = passed
        if not passed:
            all_pass = False

    result["valid"] = all_pass
    return result


def validate_formulas(formulas: list, rules: list) -> list:
    """Validate a list of formulas.

    Parameters
    ----------
    formulas:
        List of molecular formula strings.
    rules:
        List of rule names.

    Returns
    -------
    list of validation result dicts.
    """
    return [validate_formula(f, rules) for f in formulas]


@click.command()
@click.option("--input", "input_file", required=True, help="TSV file with a 'formula' column.")
@click.option("--rules", default="all",
              help="Comma-separated rules: rdbe,hc,nc,oc,sc,pc or 'all' (default: all).")
@click.option("--output", required=True, help="Output TSV with validation results.")
def main(input_file, rules, output) -> None:
    """CLI entry point."""
    rules = [r.strip() for r in rules.split(",")]

    formulas = []
    with open(input_file, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            formulas.append(row["formula"])

    results = validate_formulas(formulas, rules)

    fieldnames = list(results[0].keys()) if results else []
    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    n_valid = sum(1 for r in results if r["valid"])
    print(f"Validated {len(results)} formulas: {n_valid} passed, {len(results) - n_valid} failed.")


if __name__ == "__main__":
    main()
