"""
Intensity Distribution Reporter
================================
Report per-sample intensity statistics from a quantification matrix.

Computes mean, median, standard deviation, min, max, Q1, Q3, and count of
non-missing values for each sample column.

Usage
-----
    python intensity_distribution_reporter.py --input matrix.tsv --output intensity_stats.tsv
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


def compute_intensity_stats(matrix: np.ndarray, col_names: list) -> list:
    """Compute per-sample intensity statistics.

    Parameters
    ----------
    matrix:
        2D array (features x samples).
    col_names:
        Sample names.

    Returns
    -------
    list
        List of dicts with keys: sample, n_values, n_missing, mean, median, sd, min, max, q1, q3.
    """
    n_features = matrix.shape[0]
    results = []

    for col_idx, sample in enumerate(col_names):
        col_data = matrix[:, col_idx]
        valid = col_data[~np.isnan(col_data)]
        n_valid = len(valid)
        n_missing = n_features - n_valid

        if n_valid == 0:
            results.append({
                "sample": sample,
                "n_values": 0,
                "n_missing": n_features,
                "mean": float("nan"),
                "median": float("nan"),
                "sd": float("nan"),
                "min": float("nan"),
                "max": float("nan"),
                "q1": float("nan"),
                "q3": float("nan"),
            })
        else:
            results.append({
                "sample": sample,
                "n_values": n_valid,
                "n_missing": n_missing,
                "mean": float(np.mean(valid)),
                "median": float(np.median(valid)),
                "sd": float(np.std(valid, ddof=1)) if n_valid > 1 else 0.0,
                "min": float(np.min(valid)),
                "max": float(np.max(valid)),
                "q1": float(np.percentile(valid, 25)),
                "q3": float(np.percentile(valid, 75)),
            })

    return results


def main():
    parser = argparse.ArgumentParser(description="Per-sample intensity statistics.")
    parser.add_argument("--input", required=True, help="Input TSV matrix file")
    parser.add_argument("--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    row_ids, col_names, matrix = read_matrix(args.input)
    stats = compute_intensity_stats(matrix, col_names)

    fieldnames = ["sample", "n_values", "n_missing", "mean", "median", "sd", "min", "max", "q1", "q3"]
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for s in stats:
            row = {"sample": s["sample"], "n_values": s["n_values"], "n_missing": s["n_missing"]}
            for key in ["mean", "median", "sd", "min", "max", "q1", "q3"]:
                row[key] = f"{s[key]:.6f}" if not np.isnan(s[key]) else "NA"
            writer.writerow(row)

    print(f"Samples: {len(col_names)}")
    print(f"Features: {len(row_ids)}")
    print(f"Output written to {args.output}")


if __name__ == "__main__":
    main()
