"""
Quantification Normalizer
=========================
Normalize quantification matrices using median, quantile, or total intensity methods.

Usage
-----
    python quantification_normalizer.py --input matrix.tsv --method median --output normalized.tsv
    python quantification_normalizer.py --input matrix.tsv --method quantile --output normalized.tsv
"""

import argparse
import csv
import sys

try:
    import pyopenms as oms  # noqa: F401
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

import numpy as np


def read_matrix(filepath: str) -> tuple:
    """Read a TSV quantification matrix.

    Returns
    -------
    tuple
        (row_ids, col_names, data_matrix).
    """
    with open(filepath) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        col_names = header[1:]
        row_ids = []
        rows = []
        for row in reader:
            row_ids.append(row[0])
            rows.append([float(v) if v.strip() else 0.0 for v in row[1:]])
    return row_ids, col_names, np.array(rows, dtype=float)


def write_matrix(filepath: str, row_ids: list, col_names: list, matrix: np.ndarray) -> None:
    """Write a quantification matrix to TSV."""
    with open(filepath, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([""] + col_names)
        for i, row_id in enumerate(row_ids):
            writer.writerow([row_id] + [f"{v:.6f}" for v in matrix[i]])


def normalize_median(matrix: np.ndarray) -> np.ndarray:
    """Median normalization: shift each column so all columns have the same median.

    Parameters
    ----------
    matrix:
        2D numpy array (rows=features, cols=samples).

    Returns
    -------
    np.ndarray
        Normalized matrix.
    """
    col_medians = np.median(matrix, axis=0)
    global_median = np.median(col_medians)
    shifts = global_median - col_medians
    return matrix + shifts[np.newaxis, :]


def normalize_quantile(matrix: np.ndarray) -> np.ndarray:
    """Quantile normalization: force all columns to have the same distribution.

    Parameters
    ----------
    matrix:
        2D numpy array.

    Returns
    -------
    np.ndarray
        Quantile-normalized matrix.
    """
    n_rows, n_cols = matrix.shape
    sorted_indices = np.argsort(matrix, axis=0)
    sorted_matrix = np.sort(matrix, axis=0)
    row_means = np.mean(sorted_matrix, axis=1)

    result = np.empty_like(matrix)
    for col in range(n_cols):
        ranks = np.empty(n_rows, dtype=int)
        ranks[sorted_indices[:, col]] = np.arange(n_rows)
        result[:, col] = row_means[ranks]
    return result


def normalize_total_intensity(matrix: np.ndarray) -> np.ndarray:
    """Total intensity normalization: scale each column to the same total.

    Parameters
    ----------
    matrix:
        2D numpy array.

    Returns
    -------
    np.ndarray
        Normalized matrix.
    """
    col_sums = np.sum(matrix, axis=0)
    target_sum = np.mean(col_sums)
    scale_factors = target_sum / np.where(col_sums > 0, col_sums, 1.0)
    return matrix * scale_factors[np.newaxis, :]


def normalize(matrix: np.ndarray, method: str = "median") -> np.ndarray:
    """Normalize a quantification matrix.

    Parameters
    ----------
    matrix:
        2D numpy array.
    method:
        One of 'median', 'quantile', 'total_intensity'.

    Returns
    -------
    np.ndarray
        Normalized matrix.
    """
    method = method.lower()
    if method == "median":
        return normalize_median(matrix)
    elif method == "quantile":
        return normalize_quantile(matrix)
    elif method == "total_intensity":
        return normalize_total_intensity(matrix)
    else:
        raise ValueError(f"Unknown normalization method: '{method}'. Choose from: median, quantile, total_intensity")


def main():
    parser = argparse.ArgumentParser(description="Normalize quantification matrices.")
    parser.add_argument("--input", required=True, help="Input TSV matrix file")
    parser.add_argument("--method", required=True, choices=["median", "quantile", "total_intensity"],
                        help="Normalization method")
    parser.add_argument("--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    row_ids, col_names, matrix = read_matrix(args.input)
    normalized = normalize(matrix, method=args.method)
    write_matrix(args.output, row_ids, col_names, normalized)
    print(f"Method: {args.method}")
    print(f"Samples: {len(col_names)}")
    print(f"Features: {len(row_ids)}")
    print(f"Output written to {args.output}")


if __name__ == "__main__":
    main()
