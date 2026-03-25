"""
Sample Correlation Calculator
=============================
Compute Pearson or Spearman correlations between samples in a quantification matrix.

Usage
-----
    python sample_correlation_calculator.py --input matrix.tsv --method pearson --output correlations.tsv
    python sample_correlation_calculator.py --input matrix.tsv --method spearman --output correlations.tsv
"""

import argparse
import csv
import sys

try:
    import pyopenms as oms  # noqa: F401
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

import numpy as np
from scipy import stats


def read_matrix(filepath: str) -> tuple:
    """Read a TSV quantification matrix.

    Returns (row_ids, col_names, data_matrix).
    """
    with open(filepath) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        col_names = header[1:]
        row_ids = []
        rows = []
        for row in reader:
            row_ids.append(row[0])
            values = []
            for v in row[1:]:
                v = v.strip()
                if v == "" or v.upper() in ("NA", "NAN"):
                    values.append(np.nan)
                else:
                    values.append(float(v))
            rows.append(values)
    return row_ids, col_names, np.array(rows, dtype=float)


def compute_correlations(matrix: np.ndarray, col_names: list, method: str = "pearson") -> list:
    """Compute pairwise sample correlations.

    Parameters
    ----------
    matrix:
        2D array (features x samples).
    col_names:
        Sample names.
    method:
        'pearson' or 'spearman'.

    Returns
    -------
    list
        List of dicts with keys: sample_a, sample_b, correlation, pvalue.
    """
    method = method.lower()
    if method not in ("pearson", "spearman"):
        raise ValueError(f"Unknown method: '{method}'. Choose 'pearson' or 'spearman'.")

    n_samples = len(col_names)
    results = []

    for i in range(n_samples):
        for j in range(i, n_samples):
            col_i = matrix[:, i]
            col_j = matrix[:, j]
            # Use only rows where both values are non-NaN
            mask = ~np.isnan(col_i) & ~np.isnan(col_j)
            if np.sum(mask) < 3:
                corr = float("nan")
                pval = float("nan")
            else:
                if method == "pearson":
                    corr, pval = stats.pearsonr(col_i[mask], col_j[mask])
                else:
                    corr, pval = stats.spearmanr(col_i[mask], col_j[mask])

            results.append({
                "sample_a": col_names[i],
                "sample_b": col_names[j],
                "correlation": corr,
                "pvalue": pval,
            })

    return results


def correlation_matrix(matrix: np.ndarray, col_names: list, method: str = "pearson") -> np.ndarray:
    """Compute a full correlation matrix.

    Parameters
    ----------
    matrix:
        2D array (features x samples).
    col_names:
        Sample names.
    method:
        'pearson' or 'spearman'.

    Returns
    -------
    np.ndarray
        Symmetric correlation matrix.
    """
    pairs = compute_correlations(matrix, col_names, method)
    n = len(col_names)
    corr_mat = np.zeros((n, n))
    name_to_idx = {name: i for i, name in enumerate(col_names)}

    for p in pairs:
        i = name_to_idx[p["sample_a"]]
        j = name_to_idx[p["sample_b"]]
        corr_mat[i, j] = p["correlation"]
        corr_mat[j, i] = p["correlation"]

    return corr_mat


def main():
    parser = argparse.ArgumentParser(description="Compute sample correlations.")
    parser.add_argument("--input", required=True, help="Input TSV matrix file")
    parser.add_argument("--method", default="pearson", choices=["pearson", "spearman"],
                        help="Correlation method (default: pearson)")
    parser.add_argument("--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    row_ids, col_names, matrix = read_matrix(args.input)
    results = compute_correlations(matrix, col_names, args.method)

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample_a", "sample_b", "correlation", "pvalue"], delimiter="\t")
        writer.writeheader()
        for r in results:
            writer.writerow({
                "sample_a": r["sample_a"],
                "sample_b": r["sample_b"],
                "correlation": f"{r['correlation']:.6f}" if not np.isnan(r["correlation"]) else "NA",
                "pvalue": f"{r['pvalue']:.6e}" if not np.isnan(r["pvalue"]) else "NA",
            })

    print(f"Method: {args.method}")
    print(f"Samples: {len(col_names)}")
    print(f"Pairs: {len(results)}")
    print(f"Output written to {args.output}")


if __name__ == "__main__":
    main()
