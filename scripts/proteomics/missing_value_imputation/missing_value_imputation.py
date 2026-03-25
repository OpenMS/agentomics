"""
Missing Value Imputation
========================
Impute missing values in quantification matrices using various strategies.

Supported methods:
- MinDet: Replace missing values with the minimum detected value per column.
- MinProb: Replace missing values with random draws from a low-intensity
  Gaussian distribution (1st percentile).
- KNN: K-nearest-neighbor imputation using available features.

Usage
-----
    python missing_value_imputation.py --input matrix.tsv --method mindet --output imputed.tsv
    python missing_value_imputation.py --input matrix.tsv --method knn --k 5 --output imputed.tsv
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

    Parameters
    ----------
    filepath:
        Path to TSV file. First column is row IDs, first row is header.

    Returns
    -------
    tuple
        (row_ids, col_names, data_matrix) where data_matrix is a numpy array.
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
                if v.strip() == "" or v.strip().upper() == "NA" or v.strip().upper() == "NAN":
                    values.append(np.nan)
                else:
                    values.append(float(v))
            rows.append(values)
    return row_ids, col_names, np.array(rows, dtype=float)


def write_matrix(filepath: str, row_ids: list, col_names: list, matrix: np.ndarray) -> None:
    """Write a quantification matrix to TSV."""
    with open(filepath, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([""] + col_names)
        for i, row_id in enumerate(row_ids):
            writer.writerow([row_id] + [f"{v:.6f}" for v in matrix[i]])


def impute_mindet(matrix: np.ndarray) -> np.ndarray:
    """Impute missing values with minimum detected value per column.

    Parameters
    ----------
    matrix:
        2D numpy array with NaN for missing values.

    Returns
    -------
    np.ndarray
        Imputed matrix.
    """
    result = matrix.copy()
    for col in range(result.shape[1]):
        col_data = result[:, col]
        valid = col_data[~np.isnan(col_data)]
        if len(valid) > 0:
            min_val = np.min(valid)
        else:
            min_val = 0.0
        col_data[np.isnan(col_data)] = min_val
    return result


def impute_minprob(matrix: np.ndarray, q: float = 0.01, downshift: float = 1.8) -> np.ndarray:
    """Impute missing values from a low-intensity Gaussian distribution.

    For each column, draws from N(mean - downshift*sd, sd*0.3) where mean and
    sd are computed from the qth percentile of observed values.

    Parameters
    ----------
    matrix:
        2D numpy array with NaN for missing values.
    q:
        Quantile for determining the low-intensity distribution center.
    downshift:
        Number of standard deviations to shift the mean down.

    Returns
    -------
    np.ndarray
        Imputed matrix.
    """
    result = matrix.copy()
    rng = np.random.default_rng(42)
    for col in range(result.shape[1]):
        col_data = result[:, col]
        valid = col_data[~np.isnan(col_data)]
        if len(valid) == 0:
            continue
        mean_val = np.mean(valid)
        sd_val = np.std(valid) if len(valid) > 1 else mean_val * 0.1
        imp_mean = mean_val - downshift * sd_val
        imp_sd = sd_val * 0.3
        n_missing = np.sum(np.isnan(col_data))
        if n_missing > 0:
            imputed_vals = rng.normal(imp_mean, max(imp_sd, 1e-10), int(n_missing))
            col_data[np.isnan(col_data)] = imputed_vals
    return result


def impute_knn(matrix: np.ndarray, k: int = 5) -> np.ndarray:
    """Impute missing values using K-nearest-neighbor approach.

    For each row with missing values, find the k most similar rows (by
    Euclidean distance on shared observed features) and use their mean
    for imputation.

    Parameters
    ----------
    matrix:
        2D numpy array with NaN for missing values.
    k:
        Number of neighbors to use.

    Returns
    -------
    np.ndarray
        Imputed matrix.
    """
    result = matrix.copy()
    n_rows = result.shape[0]

    for i in range(n_rows):
        missing_mask = np.isnan(result[i])
        if not np.any(missing_mask):
            continue

        observed_mask = ~missing_mask
        if not np.any(observed_mask):
            # All missing: use column means
            for col in np.where(missing_mask)[0]:
                col_valid = result[:, col][~np.isnan(result[:, col])]
                result[i, col] = np.mean(col_valid) if len(col_valid) > 0 else 0.0
            continue

        # Find distances to other rows using shared observed features
        distances = []
        for j in range(n_rows):
            if i == j:
                continue
            shared = observed_mask & ~np.isnan(result[j])
            if np.sum(shared) == 0:
                continue
            dist = np.sqrt(np.sum((result[i, shared] - result[j, shared]) ** 2))
            distances.append((j, dist))

        if not distances:
            continue

        distances.sort(key=lambda x: x[1])
        neighbors = [idx for idx, _ in distances[:k]]

        for col in np.where(missing_mask)[0]:
            neighbor_vals = [result[j, col] for j in neighbors if not np.isnan(result[j, col])]
            if neighbor_vals:
                result[i, col] = np.mean(neighbor_vals)
            else:
                col_valid = result[:, col][~np.isnan(result[:, col])]
                result[i, col] = np.mean(col_valid) if len(col_valid) > 0 else 0.0

    return result


def impute(matrix: np.ndarray, method: str = "mindet", **kwargs) -> np.ndarray:
    """Impute missing values using the specified method.

    Parameters
    ----------
    matrix:
        2D numpy array with NaN for missing values.
    method:
        One of 'mindet', 'minprob', 'knn'.

    Returns
    -------
    np.ndarray
        Imputed matrix.
    """
    method = method.lower()
    if method == "mindet":
        return impute_mindet(matrix)
    elif method == "minprob":
        return impute_minprob(matrix, **kwargs)
    elif method == "knn":
        k = kwargs.get("k", 5)
        return impute_knn(matrix, k=k)
    else:
        raise ValueError(f"Unknown imputation method: '{method}'. Choose from: mindet, minprob, knn")


def main():
    parser = argparse.ArgumentParser(description="Impute missing values in quantification matrices.")
    parser.add_argument("--input", required=True, help="Input TSV matrix file")
    parser.add_argument("--method", required=True, choices=["mindet", "minprob", "knn"],
                        help="Imputation method")
    parser.add_argument("--k", type=int, default=5, help="Number of neighbors for KNN (default: 5)")
    parser.add_argument("--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    row_ids, col_names, matrix = read_matrix(args.input)
    n_missing_before = int(np.sum(np.isnan(matrix)))
    imputed = impute(matrix, method=args.method, k=args.k)
    n_missing_after = int(np.sum(np.isnan(imputed)))
    write_matrix(args.output, row_ids, col_names, imputed)

    print(f"Method: {args.method}")
    print(f"Missing values before: {n_missing_before}")
    print(f"Missing values after:  {n_missing_after}")
    print(f"Output written to {args.output}")


if __name__ == "__main__":
    main()
