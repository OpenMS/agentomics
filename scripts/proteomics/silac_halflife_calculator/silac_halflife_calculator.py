"""
SILAC Half-Life Calculator
===========================
Fit exponential decay to SILAC heavy/light (H/L) ratios for protein turnover
analysis.  For each protein, the tool fits the model::

    R(t) = R0 * exp(-k * t)

where *R(t)* is the H/L ratio at time *t*, *R0* is the initial ratio, and
*k* is the decay rate constant.  The half-life is ``ln(2) / k``.

Uses scipy.optimize.curve_fit for non-linear least-squares fitting.

Usage
-----
    python silac_halflife_calculator.py --input hl_ratios.tsv \
        --timepoints 0,6,12,24,48 --output halflives.tsv
"""

import argparse
import csv
import math
import sys
from typing import Dict, List, Optional

try:
    import pyopenms as oms  # noqa: F401
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

import numpy as np
from scipy.optimize import curve_fit


def exponential_decay(t: np.ndarray, r0: float, k: float) -> np.ndarray:
    """Exponential decay model: R(t) = R0 * exp(-k * t)."""
    return r0 * np.exp(-k * t)


def fit_halflife(
    timepoints: List[float], ratios: List[float]
) -> Optional[Dict[str, float]]:
    """Fit exponential decay to H/L ratios and compute half-life.

    Parameters
    ----------
    timepoints:
        Time values (e.g. hours).
    ratios:
        Corresponding H/L ratios.

    Returns
    -------
    dict or None
        ``r0``, ``k``, ``halflife``, ``r_squared`` if fit succeeds; None otherwise.
    """
    t = np.array(timepoints, dtype=float)
    r = np.array(ratios, dtype=float)

    # Filter out NaN/inf
    valid = np.isfinite(t) & np.isfinite(r) & (r > 0)
    t = t[valid]
    r = r[valid]

    if len(t) < 2:
        return None

    try:
        # Initial guesses
        r0_guess = float(r[0]) if r[0] > 0 else 1.0
        k_guess = 0.01
        popt, _ = curve_fit(
            exponential_decay, t, r,
            p0=[r0_guess, k_guess],
            bounds=([0, 0], [np.inf, np.inf]),
            maxfev=10000,
        )
        r0_fit, k_fit = popt

        if k_fit <= 0:
            return None

        halflife = math.log(2) / k_fit

        # R-squared
        r_pred = exponential_decay(t, r0_fit, k_fit)
        ss_res = np.sum((r - r_pred) ** 2)
        ss_tot = np.sum((r - np.mean(r)) ** 2)
        r_squared = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0

        return {
            "r0": r0_fit,
            "k": k_fit,
            "halflife": halflife,
            "r_squared": r_squared,
        }
    except (RuntimeError, ValueError):
        return None


def compute_halflives(
    proteins: Dict[str, List[float]], timepoints: List[float]
) -> List[Dict[str, object]]:
    """Compute half-lives for multiple proteins.

    Parameters
    ----------
    proteins:
        Mapping of protein ID to list of H/L ratios (one per timepoint).
    timepoints:
        Time values corresponding to ratio columns.

    Returns
    -------
    list of dict
        One entry per protein with fit results.
    """
    results: List[Dict[str, object]] = []
    for protein_id, ratios in proteins.items():
        fit = fit_halflife(timepoints, ratios)
        if fit is not None:
            results.append({
                "protein_id": protein_id,
                "r0": fit["r0"],
                "k": fit["k"],
                "halflife": fit["halflife"],
                "r_squared": fit["r_squared"],
                "status": "ok",
            })
        else:
            results.append({
                "protein_id": protein_id,
                "r0": float("nan"),
                "k": float("nan"),
                "halflife": float("nan"),
                "r_squared": float("nan"),
                "status": "fit_failed",
            })
    return results


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fit exponential decay to SILAC H/L ratios for protein turnover."
    )
    parser.add_argument(
        "--input", required=True,
        help="Input TSV with 'protein_id' and ratio columns (one per timepoint)",
    )
    parser.add_argument(
        "--timepoints", required=True,
        help="Comma-separated timepoint values (e.g. 0,6,12,24,48)",
    )
    parser.add_argument("--output", required=True, help="Output half-lives TSV")
    args = parser.parse_args()

    timepoints = [float(t) for t in args.timepoints.split(",")]

    proteins: Dict[str, List[float]] = {}
    with open(args.input, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        # Ratio columns are all columns except protein_id
        ratio_cols = [f for f in fields if f != "protein_id"]
        if len(ratio_cols) != len(timepoints):
            sys.exit(
                f"Number of ratio columns ({len(ratio_cols)}) does not match "
                f"number of timepoints ({len(timepoints)})"
            )
        for row in reader:
            pid = row.get("protein_id", "").strip()
            if not pid:
                continue
            ratios = []
            for col in ratio_cols:
                val = row.get(col, "").strip()
                try:
                    ratios.append(float(val))
                except (ValueError, TypeError):
                    ratios.append(float("nan"))
            proteins[pid] = ratios

    if not proteins:
        sys.exit("No proteins found in input.")

    results = compute_halflives(proteins, timepoints)

    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["protein_id", "r0", "k", "halflife", "r_squared", "status"])
        for r in results:
            writer.writerow([
                r["protein_id"],
                f"{r['r0']:.6f}" if not math.isnan(r["r0"]) else "NA",
                f"{r['k']:.6f}" if not math.isnan(r["k"]) else "NA",
                f"{r['halflife']:.4f}" if not math.isnan(r["halflife"]) else "NA",
                f"{r['r_squared']:.4f}" if not math.isnan(r["r_squared"]) else "NA",
                r["status"],
            ])

    ok_count = sum(1 for r in results if r["status"] == "ok")
    print(f"Fitted {ok_count}/{len(results)} proteins -> {args.output}")


if __name__ == "__main__":
    main()
