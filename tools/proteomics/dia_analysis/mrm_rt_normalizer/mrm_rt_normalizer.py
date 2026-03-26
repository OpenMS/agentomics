"""
MRM RT Normalizer
=================
Normalize retention times using MRMRTNormalizer outlier removal and linear fitting.

Features
--------
- Read iRT/measured RT pairs from TSV
- Outlier removal via iterative (Chauvenet) or RANSAC methods
- Linear model fitting on cleaned pairs
- Output cleaned pairs and model parameters

Usage
-----
    python mrm_rt_normalizer.py --input rt_pairs.tsv --output model.tsv --method iterative
    python mrm_rt_normalizer.py --input rt_pairs.tsv --output model.tsv --method ransac
"""

import csv
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def _fit_linear(pairs: list) -> tuple:
    """Fit a simple linear model y = slope * x + intercept.

    Parameters
    ----------
    pairs : list
        List of [x, y] pairs.

    Returns
    -------
    tuple
        (slope, intercept, r_squared).
    """
    n = len(pairs)
    if n < 2:
        return 0.0, 0.0, 0.0

    sum_x = sum(p[0] for p in pairs)
    sum_y = sum(p[1] for p in pairs)
    sum_xy = sum(p[0] * p[1] for p in pairs)
    sum_x2 = sum(p[0] ** 2 for p in pairs)

    denom = n * sum_x2 - sum_x ** 2
    if abs(denom) < 1e-12:
        return 0.0, sum_y / n if n > 0 else 0.0, 0.0

    slope = (n * sum_xy - sum_x * sum_y) / denom
    intercept = (sum_y - slope * sum_x) / n

    # R-squared
    mean_y = sum_y / n
    ss_tot = sum((p[1] - mean_y) ** 2 for p in pairs)
    ss_res = sum((p[1] - (slope * p[0] + intercept)) ** 2 for p in pairs)
    r_squared = 1.0 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

    return slope, intercept, r_squared


def normalize_mrm_rt(
    pairs_path: str,
    output_path: str,
    method: str = "iterative",
    rsq_limit: float = 0.5,
    coverage_limit: float = 0.3,
) -> dict:
    """Normalize MRM retention times using outlier removal and linear fitting.

    Reads RT pairs from TSV, applies outlier removal using MRMRTNormalizer,
    fits a linear model, and outputs cleaned pairs with model parameters.

    Parameters
    ----------
    pairs_path : str
        Path to input TSV file with columns: irt (or reference_rt) and measured_rt.
    output_path : str
        Path to output TSV file with cleaned pairs and model.
    method : str
        Outlier removal method: 'iterative' or 'ransac'.
    rsq_limit : float
        Minimum R-squared for convergence.
    coverage_limit : float
        Minimum coverage limit.

    Returns
    -------
    dict
        Dictionary with 'before_count', 'after_count', 'outliers_removed',
        'slope', 'intercept', 'r_squared'.
    """
    pairs = []
    with open(pairs_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            ref_col = "irt" if "irt" in row else "reference_rt"
            ref_rt = float(row[ref_col])
            meas_rt = float(row["measured_rt"])
            pairs.append([ref_rt, meas_rt])

    before_count = len(pairs)

    if method == "ransac" and len(pairs) >= 30:
        try:
            cleaned = oms.MRMRTNormalizer.removeOutliersRANSAC(
                pairs, rsq_limit, coverage_limit, 2000, 30.0, 15,
            )
        except RuntimeError:
            # Fall back to iterative if RANSAC fails
            cleaned = oms.MRMRTNormalizer.removeOutliersIterative(
                pairs, rsq_limit, coverage_limit, True, b"iter_jackknife",
            )
    else:
        try:
            cleaned = oms.MRMRTNormalizer.removeOutliersIterative(
                pairs, rsq_limit, coverage_limit, True, b"iter_jackknife",
            )
        except RuntimeError:
            # If iterative fails, keep all pairs
            cleaned = pairs

    after_count = len(cleaned)
    slope, intercept, r_squared = _fit_linear(cleaned)

    # Write output
    with open(output_path, "w", newline="") as fh:
        # Write model parameters as header comments
        fh.write(f"# slope={slope:.8f}\n")
        fh.write(f"# intercept={intercept:.8f}\n")
        fh.write(f"# r_squared={r_squared:.8f}\n")
        fh.write(f"# outliers_removed={before_count - after_count}\n")

        writer = csv.DictWriter(
            fh, fieldnames=["reference_rt", "measured_rt"], delimiter="\t"
        )
        writer.writeheader()
        for p in cleaned:
            writer.writerow({"reference_rt": round(p[0], 6), "measured_rt": round(p[1], 6)})

    return {
        "before_count": before_count,
        "after_count": after_count,
        "outliers_removed": before_count - after_count,
        "slope": slope,
        "intercept": intercept,
        "r_squared": r_squared,
    }


@click.command(help="Normalize MRM retention times with outlier removal.")
@click.option("--input", "pairs_path", required=True, help="Path to RT pairs TSV file.")
@click.option("--output", required=True, help="Path to output model TSV file.")
@click.option(
    "--method", type=click.Choice(["iterative", "ransac"]), default="iterative",
    help="Outlier removal method (default: iterative).",
)
@click.option("--rsq-limit", type=float, default=0.5, help="R-squared limit.")
@click.option("--coverage-limit", type=float, default=0.3, help="Coverage limit.")
def main(pairs_path, output, method, rsq_limit, coverage_limit):
    """CLI entry point."""
    result = normalize_mrm_rt(pairs_path, output, method, rsq_limit, coverage_limit)
    print(f"Before: {result['before_count']} pairs, After: {result['after_count']} pairs")
    print(f"Outliers removed: {result['outliers_removed']}")
    print(f"Model: y = {result['slope']:.6f} * x + {result['intercept']:.6f}")
    print(f"R-squared: {result['r_squared']:.6f}")
    print(f"Output: {output}")


if __name__ == "__main__":
    main()
