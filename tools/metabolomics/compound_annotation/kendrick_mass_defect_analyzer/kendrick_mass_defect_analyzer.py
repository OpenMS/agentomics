"""
Kendrick Mass Defect Analyzer
==============================
Compute Kendrick Mass Defect (KMD) for configurable base units (CH2, CF2, C2H4O)
and group features into homologous series.

Usage
-----
    python kendrick_mass_defect_analyzer.py --input features.tsv --base CH2 --output kmd.tsv
"""

import csv
import sys

import click
import pyopenms as oms

# Nominal masses for common base units
BASE_NOMINAL_MASSES = {
    "CH2": 14,
    "CF2": 50,
    "C2H4O": 44,
}


def get_base_exact_mass(base_formula: str) -> float:
    """Get the exact monoisotopic mass of a base unit formula.

    Parameters
    ----------
    base_formula:
        Molecular formula of the base unit, e.g. ``"CH2"``.

    Returns
    -------
    float: Exact monoisotopic mass.
    """
    ef = oms.EmpiricalFormula(base_formula)
    return ef.getMonoWeight()


def compute_kmd(exact_mass: float, base_formula: str) -> dict:
    """Compute Kendrick mass and Kendrick mass defect for a given exact mass.

    Parameters
    ----------
    exact_mass:
        Exact monoisotopic mass of the compound.
    base_formula:
        Base unit formula (e.g. ``"CH2"``).

    Returns
    -------
    dict with keys: exact_mass, kendrick_mass, nominal_kendrick_mass, kmd
    """
    base_exact = get_base_exact_mass(base_formula)
    nominal = BASE_NOMINAL_MASSES.get(base_formula)
    if nominal is None:
        ef = oms.EmpiricalFormula(base_formula)
        nominal = round(ef.getMonoWeight())

    kendrick_factor = nominal / base_exact
    kendrick_mass = exact_mass * kendrick_factor
    nominal_kendrick_mass = round(kendrick_mass)
    kmd = nominal_kendrick_mass - kendrick_mass

    return {
        "exact_mass": round(exact_mass, 6),
        "kendrick_mass": round(kendrick_mass, 6),
        "nominal_kendrick_mass": nominal_kendrick_mass,
        "kmd": round(kmd, 6),
    }


def compute_kmd_from_formula(formula: str, base_formula: str) -> dict:
    """Compute KMD from a molecular formula string.

    Parameters
    ----------
    formula:
        Molecular formula of the compound.
    base_formula:
        Base unit formula.

    Returns
    -------
    dict with keys: formula, exact_mass, kendrick_mass, nominal_kendrick_mass, kmd
    """
    ef = oms.EmpiricalFormula(formula)
    exact_mass = ef.getMonoWeight()
    result = compute_kmd(exact_mass, base_formula)
    result["formula"] = formula
    return result


def group_homologous_series(kmd_results: list, kmd_tolerance: float = 0.005) -> list:
    """Group features into homologous series based on similar KMD values.

    Parameters
    ----------
    kmd_results:
        List of dicts from compute_kmd / compute_kmd_from_formula.
    kmd_tolerance:
        Maximum KMD difference to assign features to the same series.

    Returns
    -------
    list of dicts with an added 'series' integer key.
    """
    if not kmd_results:
        return []

    sorted_results = sorted(kmd_results, key=lambda x: x["kmd"])
    series_id = 0
    sorted_results[0]["series"] = series_id

    for i in range(1, len(sorted_results)):
        if abs(sorted_results[i]["kmd"] - sorted_results[i - 1]["kmd"]) > kmd_tolerance:
            series_id += 1
        sorted_results[i]["series"] = series_id

    return sorted_results


@click.command()
@click.option("--input", "input_file", required=True, help="TSV file with 'mz' or 'formula' column.")
@click.option("--base", default="CH2", help="Base unit formula (default: CH2).")
@click.option("--kmd-tolerance", type=float, default=0.005,
              help="KMD tolerance for grouping homologous series (default: 0.005).")
@click.option("--output", required=True, help="Output TSV file with KMD values.")
def main(input_file, base, kmd_tolerance, output) -> None:
    """CLI entry point."""
    results = []
    with open(input_file, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        headers = reader.fieldnames or []
        for row in reader:
            if "formula" in headers:
                result = compute_kmd_from_formula(row["formula"], base)
            elif "mz" in headers:
                result = compute_kmd(float(row["mz"]), base)
            else:
                sys.exit("Input TSV must have a 'formula' or 'mz' column.")
            results.append(result)

    results = group_homologous_series(results, kmd_tolerance=kmd_tolerance)

    fieldnames = list(results[0].keys()) if results else []
    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    print(f"Wrote {len(results)} entries to {output}")


if __name__ == "__main__":
    main()
