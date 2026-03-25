"""
Biomarker Panel ROC
====================
Compute ROC curves and AUC values for individual protein biomarkers and
simple multi-marker panels.  For each protein, the tool computes a
receiver-operating-characteristic curve (sensitivity vs 1-specificity)
and the area under the curve (AUC) to evaluate discriminatory power
between case and control groups.

For multi-marker panels, a simple sum-score is used to combine markers.

Uses numpy and scipy for statistical computations.

Usage
-----
    python biomarker_panel_roc.py --input protein_quant.tsv \
        --groups case,control --output roc.tsv
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


def compute_roc(
    scores: List[float], labels: List[int]
) -> Tuple[List[float], List[float], float]:
    """Compute ROC curve and AUC for binary classification.

    Parameters
    ----------
    scores:
        Numeric scores (higher = more likely positive).
    labels:
        Binary labels (1 = positive/case, 0 = negative/control).

    Returns
    -------
    tuple
        (fpr_list, tpr_list, auc) where fpr_list and tpr_list define
        the ROC curve points and auc is the area under the curve.
    """
    scores_arr = np.array(scores, dtype=float)
    labels_arr = np.array(labels, dtype=int)

    # Sort by score descending
    order = np.argsort(-scores_arr)
    sorted_labels = labels_arr[order]

    n_pos = np.sum(labels_arr == 1)
    n_neg = np.sum(labels_arr == 0)

    if n_pos == 0 or n_neg == 0:
        return [0.0, 1.0], [0.0, 1.0], 0.5

    tp = 0
    fp = 0
    fpr_list = [0.0]
    tpr_list = [0.0]

    for label in sorted_labels:
        if label == 1:
            tp += 1
        else:
            fp += 1
        fpr_list.append(fp / n_neg)
        tpr_list.append(tp / n_pos)

    # AUC via trapezoidal rule
    auc = 0.0
    for i in range(1, len(fpr_list)):
        auc += (fpr_list[i] - fpr_list[i - 1]) * (tpr_list[i] + tpr_list[i - 1]) / 2.0

    return fpr_list, tpr_list, auc


def analyze_biomarkers(
    quant_data: Dict[str, Dict[str, float]],
    sample_groups: Dict[str, int],
) -> List[Dict[str, object]]:
    """Compute ROC/AUC for each protein.

    Parameters
    ----------
    quant_data:
        Mapping of protein_id to {sample_id: abundance}.
    sample_groups:
        Mapping of sample_id to label (1=case, 0=control).

    Returns
    -------
    list of dict
        One entry per protein with ``protein_id``, ``auc``, ``direction``.
    """
    results: List[Dict[str, object]] = []

    for protein_id, abundances in quant_data.items():
        scores = []
        labels = []
        for sample_id, label in sample_groups.items():
            if sample_id in abundances:
                scores.append(abundances[sample_id])
                labels.append(label)

        if len(scores) < 4:
            continue

        _, _, auc = compute_roc(scores, labels)

        # If AUC < 0.5, flip direction (lower values = case)
        direction = "up"
        if auc < 0.5:
            flipped_scores = [-s for s in scores]
            _, _, auc_flipped = compute_roc(flipped_scores, labels)
            auc = auc_flipped
            direction = "down"

        results.append({
            "protein_id": protein_id,
            "auc": auc,
            "direction": direction,
            "n_case": sum(1 for lab in labels if lab == 1),
            "n_control": sum(1 for lab in labels if lab == 0),
        })

    results.sort(key=lambda x: x["auc"], reverse=True)
    return results


def panel_score(
    quant_data: Dict[str, Dict[str, float]],
    sample_groups: Dict[str, int],
    top_markers: List[str],
) -> Tuple[List[float], List[float], float]:
    """Compute a combined panel ROC using sum-score of top markers.

    Parameters
    ----------
    quant_data:
        Per-protein abundances.
    sample_groups:
        Sample labels.
    top_markers:
        List of protein IDs to include in the panel.

    Returns
    -------
    tuple
        (fpr, tpr, auc) for the combined panel.
    """
    sample_scores: Dict[str, float] = {}
    for sample_id in sample_groups:
        total = 0.0
        count = 0
        for prot in top_markers:
            if prot in quant_data and sample_id in quant_data[prot]:
                total += quant_data[prot][sample_id]
                count += 1
        if count > 0:
            sample_scores[sample_id] = total / count

    scores = []
    labels = []
    for sample_id, score in sample_scores.items():
        scores.append(score)
        labels.append(sample_groups[sample_id])

    if len(scores) < 4:
        return [0.0, 1.0], [0.0, 1.0], 0.5

    return compute_roc(scores, labels)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute ROC/AUC for protein biomarker panels."
    )
    parser.add_argument(
        "--input", required=True,
        help="Input TSV: rows=proteins, columns=samples, first column=protein_id",
    )
    parser.add_argument(
        "--groups", required=True,
        help="Comma-separated group names: positive,negative (e.g. case,control)",
    )
    parser.add_argument(
        "--group-file", default=None,
        help="Optional TSV mapping sample_id to group. If not provided, "
             "column names must contain group labels.",
    )
    parser.add_argument("--output", required=True, help="Output ROC/AUC TSV")
    args = parser.parse_args()

    group_names = [g.strip() for g in args.groups.split(",")]
    if len(group_names) != 2:
        sys.exit("--groups must specify exactly two groups: positive,negative")
    pos_group, neg_group = group_names

    # Read quantification matrix
    quant_data: Dict[str, Dict[str, float]] = {}
    sample_ids: List[str] = []

    with open(args.input, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        sample_ids = [f for f in fields if f != "protein_id"]
        for row in reader:
            pid = row.get("protein_id", "").strip()
            if not pid:
                continue
            abundances: Dict[str, float] = {}
            for sid in sample_ids:
                val = row.get(sid, "").strip()
                try:
                    abundances[sid] = float(val)
                except (ValueError, TypeError):
                    pass
            quant_data[pid] = abundances

    # Determine sample groups
    sample_groups: Dict[str, int] = {}
    if args.group_file:
        with open(args.group_file, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                sid = row.get("sample_id", "").strip()
                grp = row.get("group", "").strip()
                if sid and grp == pos_group:
                    sample_groups[sid] = 1
                elif sid and grp == neg_group:
                    sample_groups[sid] = 0
    else:
        # Infer from column names containing group labels
        for sid in sample_ids:
            if pos_group in sid:
                sample_groups[sid] = 1
            elif neg_group in sid:
                sample_groups[sid] = 0

    if not sample_groups:
        sys.exit("Could not assign any samples to groups. Check --groups or --group-file.")

    results = analyze_biomarkers(quant_data, sample_groups)

    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["protein_id", "auc", "direction", "n_case", "n_control"])
        for r in results:
            writer.writerow([
                r["protein_id"], f"{r['auc']:.4f}",
                r["direction"], r["n_case"], r["n_control"],
            ])

    if results:
        print(f"Top marker: {results[0]['protein_id']} AUC={results[0]['auc']:.4f}")
    print(f"Analyzed {len(results)} proteins -> {args.output}")


if __name__ == "__main__":
    main()
