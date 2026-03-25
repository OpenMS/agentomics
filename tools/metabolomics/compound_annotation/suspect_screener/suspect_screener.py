"""
Suspect Screener
================
Match detected features against a suspect screening list by exact mass.
Features are matched within a user-defined ppm tolerance, and results are
ranked by mass error.

Usage
-----
    python suspect_screener.py --input features.tsv --suspects suspect_list.csv --ppm 5 --output matches.tsv
"""

import csv
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def exact_mass_from_formula(formula: str) -> float:
    """Compute the monoisotopic mass for a molecular formula using pyopenms.

    Parameters
    ----------
    formula:
        Empirical formula string, e.g. ``"C6H12O6"``.

    Returns
    -------
    float
        Monoisotopic mass in Da.
    """
    ef = oms.EmpiricalFormula(formula)
    return ef.getMonoWeight()


def ppm_error(observed: float, theoretical: float) -> float:
    """Calculate the mass error in ppm.

    Parameters
    ----------
    observed:
        Observed mass in Da.
    theoretical:
        Theoretical exact mass in Da.

    Returns
    -------
    float
        Signed mass error in ppm.
    """
    if theoretical == 0.0:
        return float("inf")
    return (observed - theoretical) / theoretical * 1e6


def load_features(path: str) -> list[dict]:
    """Load feature table from a TSV file.

    Expected columns: feature_id, mz, rt (retention time), intensity.
    The mz column is used for matching.

    Parameters
    ----------
    path:
        Path to TSV file.

    Returns
    -------
    list of dict
        Each dict has keys from the TSV header with numeric values parsed.
    """
    features = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            parsed = {}
            for key, val in row.items():
                try:
                    parsed[key] = float(val)
                except (ValueError, TypeError):
                    parsed[key] = val
            features.append(parsed)
    return features


def load_suspects(path: str) -> list[dict]:
    """Load suspect list from a CSV file.

    Expected columns: name, formula, exact_mass.
    If exact_mass is missing or empty, it will be computed from formula.

    Parameters
    ----------
    path:
        Path to CSV file.

    Returns
    -------
    list of dict
        Each dict has keys: name, formula, exact_mass.
    """
    suspects = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            name = row.get("name", "").strip()
            formula = row.get("formula", "").strip()
            mass_str = row.get("exact_mass", "").strip()
            if mass_str:
                exact_mass = float(mass_str)
            elif formula:
                exact_mass = exact_mass_from_formula(formula)
            else:
                continue
            suspects.append({
                "name": name,
                "formula": formula,
                "exact_mass": exact_mass,
            })
    return suspects


def screen_suspects(
    features: list[dict],
    suspects: list[dict],
    ppm_tolerance: float = 5.0,
    mz_column: str = "mz",
) -> list[dict]:
    """Match features against suspects within ppm tolerance.

    Parameters
    ----------
    features:
        List of feature dicts (must contain *mz_column*).
    suspects:
        List of suspect dicts with 'name', 'formula', 'exact_mass'.
    ppm_tolerance:
        Maximum absolute ppm error for a match.
    mz_column:
        Name of the m/z column in the feature table.

    Returns
    -------
    list of dict
        Matched results sorted by absolute ppm error, each containing
        feature info, suspect info, and the computed ppm error.
    """
    matches = []
    for feat in features:
        obs_mz = float(feat[mz_column])
        for suspect in suspects:
            error = ppm_error(obs_mz, suspect["exact_mass"])
            if abs(error) <= ppm_tolerance:
                match = {
                    "feature_id": feat.get("feature_id", ""),
                    "observed_mz": obs_mz,
                    "rt": feat.get("rt", ""),
                    "intensity": feat.get("intensity", ""),
                    "suspect_name": suspect["name"],
                    "suspect_formula": suspect["formula"],
                    "suspect_exact_mass": suspect["exact_mass"],
                    "ppm_error": round(error, 4),
                    "abs_ppm_error": round(abs(error), 4),
                }
                matches.append(match)
    matches.sort(key=lambda m: m["abs_ppm_error"])
    return matches


def write_matches(matches: list[dict], path: str) -> None:
    """Write match results to a TSV file.

    Parameters
    ----------
    matches:
        List of match dicts from :func:`screen_suspects`.
    path:
        Output TSV path.
    """
    if not matches:
        with open(path, "w") as fh:
            fh.write("# No matches found\n")
        return
    fieldnames = list(matches[0].keys())
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(matches)


@click.command()
@click.option("--input", "input_file", required=True, help="Feature table (TSV) with mz column")
@click.option("--suspects", required=True, help="Suspect list (CSV) with name, formula, exact_mass")
@click.option("--ppm", type=float, default=5.0, help="PPM tolerance (default: 5)")
@click.option("--output", required=True, help="Output matches (TSV)")
@click.option("--mz-column", default="mz", help="Name of m/z column in features (default: mz)")
def main(input_file, suspects, ppm, output, mz_column) -> None:
    """CLI entry point."""
    features = load_features(input_file)
    suspects_data = load_suspects(suspects)
    matches = screen_suspects(features, suspects_data, ppm_tolerance=ppm, mz_column=mz_column)
    write_matches(matches, output)
    print(f"Found {len(matches)} matches, written to {output}")


if __name__ == "__main__":
    main()
