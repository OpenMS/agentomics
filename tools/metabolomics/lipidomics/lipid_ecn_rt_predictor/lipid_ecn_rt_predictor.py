"""
Lipid ECN-RT Predictor
======================
Predict lipid retention times from Equivalent Carbon Number (ECN).
ECN = total_carbons - 2 * double_bonds.  A linear regression of RT vs ECN
is built per lipid class from calibration standards, then applied to predict
RTs for unknown lipids.

Usage
-----
    python lipid_ecn_rt_predictor.py --input lipids.tsv --calibration standards.tsv --output predictions.tsv
"""

import csv
import sys

import click

try:
    import pyopenms as oms  # noqa: F401
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

import numpy as np
from scipy import stats


def compute_ecn(total_carbons: int, double_bonds: int) -> int:
    """Calculate the Equivalent Carbon Number.

    Parameters
    ----------
    total_carbons:
        Total number of acyl-chain carbons.
    double_bonds:
        Total number of double bonds.

    Returns
    -------
    int
        ECN = total_carbons - 2 * double_bonds.
    """
    return total_carbons - 2 * double_bonds


def load_tsv(path: str) -> list[dict]:
    """Load a TSV file into a list of dicts.

    Expected columns vary by usage but typically include:
    lipid_class, total_carbons, double_bonds, rt.

    Parameters
    ----------
    path:
        Path to TSV file.

    Returns
    -------
    list of dict
    """
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            parsed = {}
            for key, val in row.items():
                try:
                    parsed[key] = float(val)
                except (ValueError, TypeError):
                    parsed[key] = val
            rows.append(parsed)
    return rows


def build_calibration_models(
    standards: list[dict],
) -> dict[str, dict]:
    """Build linear regression models (RT vs ECN) per lipid class.

    Parameters
    ----------
    standards:
        List of dicts with keys: lipid_class, total_carbons, double_bonds, rt.

    Returns
    -------
    dict
        Mapping from lipid_class to dict with keys: slope, intercept, r_value, std_err.
    """
    by_class: dict[str, list[tuple[int, float]]] = {}
    for row in standards:
        cls = str(row["lipid_class"])
        carbons = int(row["total_carbons"])
        db = int(row["double_bonds"])
        rt = float(row["rt"])
        ecn = compute_ecn(carbons, db)
        by_class.setdefault(cls, []).append((ecn, rt))

    models = {}
    for cls, points in by_class.items():
        if len(points) < 2:
            continue
        ecns = np.array([p[0] for p in points], dtype=float)
        rts = np.array([p[1] for p in points], dtype=float)
        result = stats.linregress(ecns, rts)
        models[cls] = {
            "slope": result.slope,
            "intercept": result.intercept,
            "r_value": result.rvalue,
            "std_err": result.stderr,
        }
    return models


def predict_rt(
    lipids: list[dict],
    models: dict[str, dict],
) -> list[dict]:
    """Predict RT for lipids using calibration models.

    Parameters
    ----------
    lipids:
        List of dicts with keys: lipid_class, total_carbons, double_bonds.
    models:
        Calibration models from :func:`build_calibration_models`.

    Returns
    -------
    list of dict
        Input rows augmented with ecn, predicted_rt, and model_r_value.
    """
    results = []
    for row in lipids:
        cls = str(row["lipid_class"])
        carbons = int(row["total_carbons"])
        db = int(row["double_bonds"])
        ecn = compute_ecn(carbons, db)
        result = dict(row)
        result["ecn"] = ecn
        if cls in models:
            model = models[cls]
            predicted = model["slope"] * ecn + model["intercept"]
            result["predicted_rt"] = round(predicted, 4)
            result["model_r_value"] = round(model["r_value"], 4)
        else:
            result["predicted_rt"] = "N/A"
            result["model_r_value"] = "N/A"
        results.append(result)
    return results


def write_predictions(predictions: list[dict], path: str) -> None:
    """Write prediction results to a TSV file.

    Parameters
    ----------
    predictions:
        List of prediction dicts.
    path:
        Output TSV path.
    """
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
@click.option("--input", "input_file", required=True,
              help="Lipid table (TSV) with lipid_class, total_carbons, double_bonds")
@click.option("--calibration", required=True,
              help="Standards table (TSV) with lipid_class, total_carbons, double_bonds, rt")
@click.option("--output", required=True, help="Output predictions (TSV)")
def main(input_file, calibration, output) -> None:
    """CLI entry point."""
    lipids = load_tsv(input_file)
    standards = load_tsv(calibration)
    models = build_calibration_models(standards)
    predictions = predict_rt(lipids, models)
    write_predictions(predictions, output)
    print(f"Predicted RT for {len(predictions)} lipids, written to {output}")


if __name__ == "__main__":
    main()
