"""
Differential Expression Tester
==============================
Perform t-tests with Benjamini-Hochberg correction on quantification matrices.

Reads a quantification matrix and an experimental design file that maps samples
to conditions, then computes per-feature differential expression statistics.

Usage
-----
    python differential_expression_tester.py --input matrix.tsv --design design.tsv --test ttest --output de.tsv
    python differential_expression_tester.py --input matrix.tsv --design design.tsv --test welch --output de.tsv
"""

import argparse
import csv
import math
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


def read_design(filepath: str) -> dict:
    """Read experimental design file mapping samples to conditions.

    Expected format: TSV with columns 'sample' and 'condition'.

    Returns
    -------
    dict
        {sample_name: condition_name}
    """
    design = {}
    with open(filepath) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            design[row["sample"]] = row["condition"]
    return design


def benjamini_hochberg(pvalues: list) -> list:
    """Apply Benjamini-Hochberg FDR correction.

    Parameters
    ----------
    pvalues:
        List of p-values (may contain NaN).

    Returns
    -------
    list
        Adjusted p-values.
    """
    n = len(pvalues)
    valid_indices = [i for i in range(n) if not math.isnan(pvalues[i])]
    if not valid_indices:
        return [float("nan")] * n

    sorted_valid = sorted(valid_indices, key=lambda i: pvalues[i])
    adjusted = [float("nan")] * n
    m = len(sorted_valid)

    prev = 1.0
    for rank_idx in range(m - 1, -1, -1):
        i = sorted_valid[rank_idx]
        rank = rank_idx + 1
        adj = min(pvalues[i] * m / rank, prev)
        adj = min(adj, 1.0)
        adjusted[i] = adj
        prev = adj

    return adjusted


def differential_expression(
    matrix: np.ndarray, row_ids: list, col_names: list, design: dict, test: str = "ttest"
) -> list:
    """Compute differential expression statistics.

    Parameters
    ----------
    matrix:
        2D array (features x samples).
    row_ids:
        Feature identifiers.
    col_names:
        Sample names matching matrix columns.
    design:
        {sample: condition} mapping. Expects exactly two conditions.
    test:
        Statistical test: 'ttest' (equal variance) or 'welch' (unequal variance).

    Returns
    -------
    list
        List of dicts with keys: feature, log2fc, pvalue, adj_pvalue.
    """
    conditions = sorted(set(design.values()))
    if len(conditions) != 2:
        raise ValueError(f"Exactly 2 conditions required, got {len(conditions)}: {conditions}")

    cond_a, cond_b = conditions
    idx_a = [i for i, s in enumerate(col_names) if design.get(s) == cond_a]
    idx_b = [i for i, s in enumerate(col_names) if design.get(s) == cond_b]

    if not idx_a or not idx_b:
        raise ValueError("No samples found for one or both conditions.")

    equal_var = test.lower() == "ttest"

    results = []
    pvalues = []
    for row_idx in range(len(row_ids)):
        vals_a = matrix[row_idx, idx_a]
        vals_b = matrix[row_idx, idx_b]
        valid_a = vals_a[~np.isnan(vals_a)]
        valid_b = vals_b[~np.isnan(vals_b)]

        if len(valid_a) < 2 or len(valid_b) < 2:
            log2fc = float("nan")
            pval = float("nan")
        else:
            mean_a = np.mean(valid_a)
            mean_b = np.mean(valid_b)
            if mean_a > 0 and mean_b > 0:
                log2fc = np.log2(mean_b / mean_a)
            else:
                log2fc = float("nan")
            _, pval = stats.ttest_ind(valid_a, valid_b, equal_var=equal_var)

        results.append({
            "feature": row_ids[row_idx],
            "log2fc": log2fc,
            "pvalue": pval,
        })
        pvalues.append(pval)

    adj_pvalues = benjamini_hochberg(pvalues)
    for i, r in enumerate(results):
        r["adj_pvalue"] = adj_pvalues[i]

    return results


def main():
    parser = argparse.ArgumentParser(description="T-test + BH correction on quantification matrices.")
    parser.add_argument("--input", required=True, help="Input TSV matrix file")
    parser.add_argument("--design", required=True, help="Experimental design TSV (columns: sample, condition)")
    parser.add_argument("--test", default="ttest", choices=["ttest", "welch"], help="Test type (default: ttest)")
    parser.add_argument("--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    row_ids, col_names, matrix = read_matrix(args.input)
    design = read_design(args.design)
    results = differential_expression(matrix, row_ids, col_names, design, test=args.test)

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["feature", "log2fc", "pvalue", "adj_pvalue"], delimiter="\t")
        writer.writeheader()
        for r in results:
            writer.writerow({
                "feature": r["feature"],
                "log2fc": f"{r['log2fc']:.6f}" if not math.isnan(r["log2fc"]) else "NA",
                "pvalue": f"{r['pvalue']:.6e}" if not math.isnan(r["pvalue"]) else "NA",
                "adj_pvalue": f"{r['adj_pvalue']:.6e}" if not math.isnan(r["adj_pvalue"]) else "NA",
            })

    n_sig = sum(1 for r in results if not math.isnan(r["adj_pvalue"]) and r["adj_pvalue"] < 0.05)
    print(f"Test: {args.test}")
    print(f"Features tested: {len(results)}")
    print(f"Significant (adj_pvalue < 0.05): {n_sig}")
    print(f"Output written to {args.output}")


if __name__ == "__main__":
    main()
