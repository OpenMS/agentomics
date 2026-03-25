"""
iRT Calculator
===============
Convert observed retention times to indexed retention times (iRT) using reference peptides.

Features
--------
- Linear regression fitting between observed RT and known iRT values
- Convert observed RTs to iRT scale
- Report R-squared and regression parameters
- Support for custom reference peptide sets

Usage
-----
    python irt_calculator.py --input identifications.tsv --reference irt_standards.tsv --output irt_converted.tsv
"""

import argparse
import csv
import json
import sys

try:
    import pyopenms as oms  # noqa: F401
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def linear_regression(x: list, y: list) -> tuple:
    """Fit a simple linear regression y = slope * x + intercept.

    Parameters
    ----------
    x : list
        Independent variable values.
    y : list
        Dependent variable values.

    Returns
    -------
    tuple
        (slope, intercept, r_squared).
    """
    n = len(x)
    if n < 2:
        return 0.0, 0.0, 0.0

    sum_x = sum(x)
    sum_y = sum(y)
    sum_xy = sum(xi * yi for xi, yi in zip(x, y))
    sum_x2 = sum(xi * xi for xi in x)

    denom = n * sum_x2 - sum_x * sum_x
    if denom == 0:
        return 0.0, 0.0, 0.0

    slope = (n * sum_xy - sum_x * sum_y) / denom
    intercept = (sum_y - slope * sum_x) / n

    # R-squared
    ss_res = sum((yi - (slope * xi + intercept)) ** 2 for xi, yi in zip(x, y))
    mean_y = sum_y / n
    ss_tot = sum((yi - mean_y) ** 2 for yi in y)
    r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

    return round(slope, 8), round(intercept, 8), round(r_squared, 6)


def fit_irt_model(reference_data: list) -> dict:
    """Fit an iRT conversion model from reference peptide data.

    Parameters
    ----------
    reference_data : list
        List of dicts with 'sequence', 'observed_rt', and 'irt' keys.

    Returns
    -------
    dict
        Model parameters including slope, intercept, r_squared.
    """
    observed = [float(r["observed_rt"]) for r in reference_data]
    irt_values = [float(r["irt"]) for r in reference_data]

    slope, intercept, r_squared = linear_regression(observed, irt_values)

    return {
        "slope": slope,
        "intercept": intercept,
        "r_squared": r_squared,
        "n_reference_peptides": len(reference_data),
    }


def convert_to_irt(observed_rt: float, slope: float, intercept: float) -> float:
    """Convert an observed RT to iRT using the fitted model.

    Parameters
    ----------
    observed_rt : float
        Observed retention time.
    slope : float
        Regression slope.
    intercept : float
        Regression intercept.

    Returns
    -------
    float
        Predicted iRT value.
    """
    return round(slope * observed_rt + intercept, 4)


def process_identifications(identifications: list, model: dict) -> list:
    """Convert observed RTs to iRT for a list of identifications.

    Parameters
    ----------
    identifications : list
        List of dicts with 'sequence' and 'rt' keys.
    model : dict
        Fitted model with 'slope' and 'intercept'.

    Returns
    -------
    list
        List of dicts with added 'irt' key.
    """
    results = []
    for ident in identifications:
        rt = float(ident.get("rt", 0))
        irt = convert_to_irt(rt, model["slope"], model["intercept"])
        results.append({
            "sequence": ident.get("sequence", ""),
            "observed_rt": rt,
            "irt": irt,
        })
    return results


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(description="Convert observed RT to indexed RT (iRT).")
    parser.add_argument("--input", required=True, help="TSV with sequence and rt columns.")
    parser.add_argument("--reference", required=True, help="TSV with sequence, observed_rt, irt columns.")
    parser.add_argument("--output", help="Output file (.tsv or .json).")
    args = parser.parse_args()

    # Load reference peptides
    reference_data = []
    with open(args.reference) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            reference_data.append(row)

    model = fit_irt_model(reference_data)
    print(f"Model: iRT = {model['slope']} * RT + {model['intercept']} (R2={model['r_squared']})")

    # Load identifications
    identifications = []
    with open(args.input) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            identifications.append(row)

    results = process_identifications(identifications, model)

    if args.output:
        if args.output.endswith(".json"):
            with open(args.output, "w") as fh:
                json.dump({"model": model, "results": results}, fh, indent=2)
        else:
            with open(args.output, "w", newline="") as fh:
                writer = csv.DictWriter(fh, fieldnames=["sequence", "observed_rt", "irt"], delimiter="\t")
                writer.writeheader()
                writer.writerows(results)
        print(f"Results written to {args.output}")
    else:
        for r in results:
            print(f"{r['sequence']}\t{r['observed_rt']}\t{r['irt']}")


if __name__ == "__main__":
    main()
