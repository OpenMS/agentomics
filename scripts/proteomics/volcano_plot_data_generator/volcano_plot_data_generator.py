"""
Volcano Plot Data Generator
============================
Generate volcano plot data from differential expression results.

Annotates features as 'up', 'down', or 'ns' (not significant) based on
fold-change and p-value thresholds, and computes -log10(p-value) for plotting.

Usage
-----
    python volcano_plot_data_generator.py --input de_results.tsv --fc-threshold 1.0 --pvalue 0.05 --output volcano.tsv
"""

import argparse
import csv
import math
import sys

try:
    import pyopenms as oms  # noqa: F401
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def read_de_results(filepath: str) -> list:
    """Read differential expression results.

    Expected columns: feature, log2fc, pvalue (or adj_pvalue).

    Returns
    -------
    list
        List of dicts with keys: feature, log2fc, pvalue.
    """
    results = []
    with open(filepath) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            feature = row.get("feature", row.get("protein", row.get("peptide", "")))
            log2fc_str = row.get("log2fc", "NA")
            pval_str = row.get("adj_pvalue", row.get("pvalue", "NA"))

            try:
                log2fc = float(log2fc_str)
            except (ValueError, TypeError):
                log2fc = float("nan")
            try:
                pval = float(pval_str)
            except (ValueError, TypeError):
                pval = float("nan")

            results.append({"feature": feature, "log2fc": log2fc, "pvalue": pval})
    return results


def generate_volcano_data(
    de_results: list, fc_threshold: float = 1.0, pvalue_threshold: float = 0.05
) -> list:
    """Annotate DE results for volcano plotting.

    Parameters
    ----------
    de_results:
        List of dicts with keys: feature, log2fc, pvalue.
    fc_threshold:
        Absolute log2 fold-change threshold.
    pvalue_threshold:
        P-value significance threshold.

    Returns
    -------
    list
        List of dicts with keys: feature, log2fc, pvalue, neg_log10_pvalue, regulation.
    """
    volcano = []
    for r in de_results:
        log2fc = r["log2fc"]
        pval = r["pvalue"]

        if math.isnan(log2fc) or math.isnan(pval):
            neg_log10_p = float("nan")
            regulation = "ns"
        else:
            neg_log10_p = -math.log10(pval) if pval > 0 else float("inf")
            if pval < pvalue_threshold and log2fc > fc_threshold:
                regulation = "up"
            elif pval < pvalue_threshold and log2fc < -fc_threshold:
                regulation = "down"
            else:
                regulation = "ns"

        volcano.append({
            "feature": r["feature"],
            "log2fc": log2fc,
            "pvalue": pval,
            "neg_log10_pvalue": neg_log10_p,
            "regulation": regulation,
        })
    return volcano


def summarize_volcano(volcano_data: list) -> dict:
    """Count features by regulation status.

    Returns
    -------
    dict
        {up: int, down: int, ns: int}
    """
    counts = {"up": 0, "down": 0, "ns": 0}
    for v in volcano_data:
        counts[v["regulation"]] += 1
    return counts


def main():
    parser = argparse.ArgumentParser(description="Generate volcano plot data from DE results.")
    parser.add_argument("--input", required=True, help="Input DE results TSV")
    parser.add_argument("--fc-threshold", type=float, default=1.0, help="Log2 fold-change threshold (default: 1.0)")
    parser.add_argument("--pvalue", type=float, default=0.05, help="P-value threshold (default: 0.05)")
    parser.add_argument("--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    de_results = read_de_results(args.input)
    volcano = generate_volcano_data(de_results, fc_threshold=args.fc_threshold, pvalue_threshold=args.pvalue)

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["feature", "log2fc", "pvalue", "neg_log10_pvalue", "regulation"],
            delimiter="\t",
        )
        writer.writeheader()
        for v in volcano:
            writer.writerow({
                "feature": v["feature"],
                "log2fc": f"{v['log2fc']:.6f}" if not math.isnan(v["log2fc"]) else "NA",
                "pvalue": f"{v['pvalue']:.6e}" if not math.isnan(v["pvalue"]) else "NA",
                "neg_log10_pvalue": (
                    f"{v['neg_log10_pvalue']:.4f}" if not math.isnan(v["neg_log10_pvalue"]) else "NA"
                ),
                "regulation": v["regulation"],
            })

    counts = summarize_volcano(volcano)
    print(f"Total features: {len(volcano)}")
    print(f"Up-regulated: {counts['up']}")
    print(f"Down-regulated: {counts['down']}")
    print(f"Not significant: {counts['ns']}")
    print(f"Output written to {args.output}")


if __name__ == "__main__":
    main()
