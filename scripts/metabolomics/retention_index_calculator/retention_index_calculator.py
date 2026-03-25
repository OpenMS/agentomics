"""
Retention Index Calculator
===========================
Calculate Kovats retention indices from alkane standard retention times.

Given a set of n-alkane standards with known carbon numbers and their
observed retention times, compute retention indices for unknown compounds.

Usage
-----
    python retention_index_calculator.py --input features.tsv --standards alkanes.tsv --output ri.tsv
"""

import argparse
import csv
import math
import sys

try:
    import pyopenms as oms  # noqa: F401
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def load_standards(path: str) -> list[tuple[int, float]]:
    """Load alkane standards from a TSV file.

    Parameters
    ----------
    path:
        TSV with columns: carbon_number, rt (retention time in seconds).

    Returns
    -------
    list[tuple[int, float]]
        Sorted list of (carbon_number, rt) tuples.
    """
    standards = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            standards.append((int(row["carbon_number"]), float(row["rt"])))
    standards.sort(key=lambda x: x[1])
    return standards


def calculate_kovats_ri(
    rt: float,
    standards: list[tuple[int, float]],
) -> float | None:
    """Calculate Kovats retention index for a given retention time.

    Parameters
    ----------
    rt:
        Retention time of the unknown compound (seconds).
    standards:
        Sorted list of (carbon_number, rt) for alkane standards.

    Returns
    -------
    float or None
        Kovats retention index, or None if RT is outside the standard range.
    """
    if len(standards) < 2:
        return None

    # Find bracketing standards
    for i in range(len(standards) - 1):
        cn_z, rt_z = standards[i]
        cn_z1, rt_z1 = standards[i + 1]

        if rt_z <= rt <= rt_z1:
            if rt_z1 == rt_z:
                return float(cn_z * 100)
            # Kovats RI (isothermal): RI = 100 * [z + (log(rt_x) - log(rt_z)) / (log(rt_z1) - log(rt_z))]
            if rt_z > 0 and rt > 0 and rt_z1 > 0:
                ri = 100.0 * (cn_z + (math.log(rt) - math.log(rt_z)) / (math.log(rt_z1) - math.log(rt_z)))
                return round(ri, 2)
            else:
                # Linear interpolation fallback
                ri = 100.0 * (cn_z + (rt - rt_z) / (rt_z1 - rt_z))
                return round(ri, 2)

    return None


def calculate_all_ri(
    features: list[dict],
    standards: list[tuple[int, float]],
) -> list[dict]:
    """Calculate retention indices for all features.

    Parameters
    ----------
    features:
        List of dicts with at least key ``rt``.
    standards:
        Alkane standard reference points.

    Returns
    -------
    list[dict]
        Each feature dict augmented with ``retention_index``.
    """
    results = []
    for feat in features:
        rt = float(feat["rt"])
        ri = calculate_kovats_ri(rt, standards)
        feat_copy = dict(feat)
        feat_copy["retention_index"] = ri if ri is not None else ""
        results.append(feat_copy)
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Calculate Kovats retention indices from alkane standards."
    )
    parser.add_argument("--input", required=True, metavar="FILE", help="Features TSV (must have rt column)")
    parser.add_argument("--standards", required=True, metavar="FILE", help="Alkane standards TSV")
    parser.add_argument("--output", required=True, metavar="FILE", help="Output TSV with retention indices")
    args = parser.parse_args()

    features = []
    with open(args.input) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            features.append(row)

    standards = load_standards(args.standards)
    results = calculate_all_ri(features, standards)

    base_fields = list(features[0].keys()) if features else ["mz", "rt", "intensity"]
    fieldnames = base_fields + ["retention_index"]
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    n_annotated = sum(1 for r in results if r["retention_index"] != "")
    print(f"RI calculated for {n_annotated}/{len(results)} features, written to {args.output}")


if __name__ == "__main__":
    main()
