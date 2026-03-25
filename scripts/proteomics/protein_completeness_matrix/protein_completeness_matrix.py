"""
Protein Completeness Matrix
=============================
Compute data completeness per protein and per sample from a quantification
matrix.  Reports the fraction of non-missing values for each protein across
samples and for each sample across proteins, and optionally filters proteins
below a minimum completeness threshold.

Usage
-----
    python protein_completeness_matrix.py --input quant_matrix.tsv \
        --min-completeness 0.5 --output completeness.tsv
"""

import argparse
import csv
import sys
from typing import Dict, List, Tuple

try:
    import pyopenms as oms  # noqa: F401
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

import numpy as np


def load_quant_matrix(input_path: str) -> Tuple[List[str], List[str], np.ndarray]:
    """Load a protein quantification matrix from TSV.

    The first column is expected to be ``protein_id`` and remaining columns
    are sample IDs.  Missing values (empty strings, ``NA``, ``NaN``) are
    treated as missing.

    Parameters
    ----------
    input_path:
        Path to input TSV.

    Returns
    -------
    tuple
        (protein_ids, sample_ids, data_matrix) where data_matrix has shape
        (n_proteins, n_samples) with NaN for missing values.
    """
    protein_ids: List[str] = []
    rows: List[List[float]] = []

    with open(input_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        sample_ids = [f for f in fields if f != "protein_id"]

        for row in reader:
            pid = row.get("protein_id", "").strip()
            if not pid:
                continue
            protein_ids.append(pid)
            values: List[float] = []
            for sid in sample_ids:
                val = row.get(sid, "").strip()
                if val in ("", "NA", "NaN", "nan", "null"):
                    values.append(float("nan"))
                else:
                    try:
                        v = float(val)
                        values.append(v if v > 0 else float("nan"))
                    except (ValueError, TypeError):
                        values.append(float("nan"))
            rows.append(values)

    data = np.array(rows, dtype=float) if rows else np.empty((0, len(sample_ids)))
    return protein_ids, sample_ids, data


def compute_protein_completeness(
    data: np.ndarray,
) -> np.ndarray:
    """Compute completeness (fraction of non-NaN values) per protein (row).

    Parameters
    ----------
    data:
        Matrix of shape (n_proteins, n_samples).

    Returns
    -------
    numpy.ndarray
        Array of shape (n_proteins,) with completeness fractions.
    """
    n_samples = data.shape[1]
    if n_samples == 0:
        return np.zeros(data.shape[0])
    non_missing = np.sum(~np.isnan(data), axis=1)
    return non_missing / n_samples


def compute_sample_completeness(
    data: np.ndarray,
) -> np.ndarray:
    """Compute completeness (fraction of non-NaN values) per sample (column).

    Parameters
    ----------
    data:
        Matrix of shape (n_proteins, n_samples).

    Returns
    -------
    numpy.ndarray
        Array of shape (n_samples,) with completeness fractions.
    """
    n_proteins = data.shape[0]
    if n_proteins == 0:
        return np.zeros(data.shape[1])
    non_missing = np.sum(~np.isnan(data), axis=0)
    return non_missing / n_proteins


def filter_by_completeness(
    protein_ids: List[str],
    data: np.ndarray,
    protein_completeness: np.ndarray,
    min_completeness: float,
) -> Tuple[List[str], np.ndarray]:
    """Filter proteins by minimum completeness threshold.

    Parameters
    ----------
    protein_ids:
        Protein identifiers.
    data:
        Data matrix.
    protein_completeness:
        Per-protein completeness values.
    min_completeness:
        Minimum fraction required to keep a protein.

    Returns
    -------
    tuple
        (filtered_ids, filtered_data)
    """
    mask = protein_completeness >= min_completeness
    filtered_ids = [pid for pid, keep in zip(protein_ids, mask) if keep]
    filtered_data = data[mask]
    return filtered_ids, filtered_data


def completeness_summary(
    protein_completeness: np.ndarray,
    sample_completeness: np.ndarray,
    data: np.ndarray,
) -> Dict[str, object]:
    """Compute overall completeness summary.

    Returns
    -------
    dict
        Summary statistics.
    """
    total_cells = data.size
    non_missing = int(np.sum(~np.isnan(data)))
    return {
        "total_proteins": data.shape[0],
        "total_samples": data.shape[1],
        "total_cells": total_cells,
        "non_missing_cells": non_missing,
        "overall_completeness": non_missing / total_cells if total_cells > 0 else 0.0,
        "mean_protein_completeness": float(np.mean(protein_completeness)) if len(protein_completeness) > 0 else 0.0,
        "median_protein_completeness": float(np.median(protein_completeness)) if len(protein_completeness) > 0 else 0.0,
        "mean_sample_completeness": float(np.mean(sample_completeness)) if len(sample_completeness) > 0 else 0.0,
        "median_sample_completeness": float(np.median(sample_completeness)) if len(sample_completeness) > 0 else 0.0,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute data completeness per protein and sample."
    )
    parser.add_argument("--input", required=True, help="Input quantification matrix TSV")
    parser.add_argument(
        "--min-completeness", type=float, default=0.0,
        help="Minimum protein completeness to retain (default: 0.0 = keep all)",
    )
    parser.add_argument("--output", required=True, help="Output completeness TSV")
    args = parser.parse_args()

    protein_ids, sample_ids, data = load_quant_matrix(args.input)
    if len(protein_ids) == 0:
        sys.exit("No proteins found in input.")

    prot_comp = compute_protein_completeness(data)
    samp_comp = compute_sample_completeness(data)
    summary = completeness_summary(prot_comp, samp_comp, data)

    # Filter
    filtered_ids, filtered_data = filter_by_completeness(
        protein_ids, data, prot_comp, args.min_completeness
    )

    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        # Per-protein completeness
        writer.writerow(["protein_id", "completeness", "n_present", "n_total"])
        for i, pid in enumerate(protein_ids):
            n_present = int(np.sum(~np.isnan(data[i])))
            writer.writerow([pid, f"{prot_comp[i]:.4f}", n_present, data.shape[1]])
        writer.writerow([])

        # Per-sample completeness
        writer.writerow(["sample_id", "completeness", "n_present", "n_total"])
        for j, sid in enumerate(sample_ids):
            n_present = int(np.sum(~np.isnan(data[:, j])))
            writer.writerow([sid, f"{samp_comp[j]:.4f}", n_present, data.shape[0]])
        writer.writerow([])

        # Summary
        writer.writerow(["metric", "value"])
        for key, val in summary.items():
            if isinstance(val, float):
                writer.writerow([key, f"{val:.4f}"])
            else:
                writer.writerow([key, val])

        if args.min_completeness > 0:
            writer.writerow([])
            writer.writerow(["filter_threshold", args.min_completeness])
            writer.writerow(["proteins_retained", len(filtered_ids)])
            writer.writerow(["proteins_removed", len(protein_ids) - len(filtered_ids)])

    print(f"Proteins: {len(protein_ids)}, overall completeness: {summary['overall_completeness']:.1%}")
    if args.min_completeness > 0:
        print(
            f"After filtering (>={args.min_completeness:.0%}): "
            f"{len(filtered_ids)}/{len(protein_ids)} proteins retained"
        )
    print(f"Output -> {args.output}")


if __name__ == "__main__":
    main()
