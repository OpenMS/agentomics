"""
Coefficient of Variation Calculator
====================================
Calculate CV% (coefficient of variation) across replicates for each feature.

Reads a quantification matrix and a groups file that assigns samples to groups,
then computes the CV within each group for each feature.

Usage
-----
    python coefficient_of_variation_calculator.py --input matrix.tsv --groups groups.tsv --output cv_report.tsv
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


def read_groups(filepath: str) -> dict:
    """Read a groups file mapping samples to groups.

    Expected format: TSV with columns 'sample' and 'group'.

    Returns
    -------
    dict
        {sample_name: group_name}
    """
    groups = {}
    with open(filepath) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            groups[row["sample"]] = row["group"]
    return groups


def calculate_cv(matrix: np.ndarray, row_ids: list, col_names: list, groups: dict) -> list:
    """Calculate CV% for each feature within each group.

    Parameters
    ----------
    matrix:
        2D array (features x samples).
    row_ids:
        Feature identifiers.
    col_names:
        Sample names.
    groups:
        {sample: group} mapping.

    Returns
    -------
    list
        List of dicts with keys: feature, group, mean, sd, cv_percent, n_values.
    """
    group_names = sorted(set(groups.values()))
    group_indices = {}
    for g in group_names:
        group_indices[g] = [i for i, s in enumerate(col_names) if groups.get(s) == g]

    results = []
    for row_idx, feature in enumerate(row_ids):
        for group in group_names:
            indices = group_indices[group]
            if not indices:
                continue
            values = matrix[row_idx, indices]
            valid = values[~np.isnan(values)]
            n_valid = len(valid)
            if n_valid < 2:
                mean_val = np.mean(valid) if n_valid > 0 else float("nan")
                results.append({
                    "feature": feature,
                    "group": group,
                    "mean": mean_val,
                    "sd": float("nan"),
                    "cv_percent": float("nan"),
                    "n_values": n_valid,
                })
            else:
                mean_val = np.mean(valid)
                sd_val = np.std(valid, ddof=1)
                cv = (sd_val / mean_val * 100.0) if mean_val != 0 else float("nan")
                results.append({
                    "feature": feature,
                    "group": group,
                    "mean": mean_val,
                    "sd": sd_val,
                    "cv_percent": cv,
                    "n_values": n_valid,
                })
    return results


def main():
    parser = argparse.ArgumentParser(description="Calculate CV% across replicates.")
    parser.add_argument("--input", required=True, help="Input TSV matrix file")
    parser.add_argument("--groups", required=True, help="Groups TSV (columns: sample, group)")
    parser.add_argument("--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    row_ids, col_names, matrix = read_matrix(args.input)
    groups = read_groups(args.groups)
    results = calculate_cv(matrix, row_ids, col_names, groups)

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["feature", "group", "mean", "sd", "cv_percent", "n_values"],
            delimiter="\t",
        )
        writer.writeheader()
        for r in results:
            writer.writerow({
                "feature": r["feature"],
                "group": r["group"],
                "mean": f"{r['mean']:.6f}" if not np.isnan(r["mean"]) else "NA",
                "sd": f"{r['sd']:.6f}" if not np.isnan(r["sd"]) else "NA",
                "cv_percent": f"{r['cv_percent']:.2f}" if not np.isnan(r["cv_percent"]) else "NA",
                "n_values": r["n_values"],
            })

    valid_cvs = [r["cv_percent"] for r in results if not np.isnan(r["cv_percent"])]
    if valid_cvs:
        print(f"Features: {len(row_ids)}")
        print(f"Groups: {len(set(groups.values()))}")
        print(f"Median CV%: {np.median(valid_cvs):.2f}")
        print(f"Mean CV%: {np.mean(valid_cvs):.2f}")
    print(f"Output written to {args.output}")


if __name__ == "__main__":
    main()
