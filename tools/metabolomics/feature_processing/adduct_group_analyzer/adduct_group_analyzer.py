"""
Adduct Group Analyzer
======================
Group features that likely originate from the same compound by detecting
adduct relationships based on m/z differences and RT proximity.

Features with matching m/z differences (within tolerance) for known
adduct pairs that also co-elute are assigned to the same group.

Usage
-----
    python adduct_group_analyzer.py --input features.tsv --rt-tolerance 5 --output groups.tsv
"""

import csv

import click
import pyopenms as oms  # noqa: F401

PROTON = 1.007276

# Known adduct mass differences relative to [M+H]+ (positive mode)
ADDUCT_DIFFS = {
    "[M+Na]+ vs [M+H]+": 22.989218 - PROTON,
    "[M+K]+ vs [M+H]+": 38.963158 - PROTON,
    "[M+NH4]+ vs [M+H]+": 18.034164 - PROTON,
    "[M+2H]2+ vs [M+H]+": None,  # charge-state relationship, handled separately
    "[M+Na]+ vs [M+NH4]+": 22.989218 - 18.034164,
    "[M+K]+ vs [M+Na]+": 38.963158 - 22.989218,
}


def find_adduct_groups(
    features: list[dict],
    rt_tolerance: float = 5.0,
    mz_tolerance_da: float = 0.01,
) -> list[dict]:
    """Group features by adduct relationships.

    Parameters
    ----------
    features:
        List of dicts with keys: feature_id, mz, rt.
    rt_tolerance:
        Maximum RT difference in seconds for co-elution.
    mz_tolerance_da:
        Tolerance for matching adduct m/z differences in Da.

    Returns
    -------
    list[dict]
        Each dict has: feature_id, mz, rt, group_id, adduct_annotation.
    """
    n = len(features)
    group_ids = list(range(n))
    annotations = ["" for _ in range(n)]

    def find_root(i: int) -> int:
        while group_ids[i] != i:
            group_ids[i] = group_ids[group_ids[i]]
            i = group_ids[i]
        return i

    def union(i: int, j: int):
        ri, rj = find_root(i), find_root(j)
        if ri != rj:
            group_ids[rj] = ri

    active_diffs = {name: diff for name, diff in ADDUCT_DIFFS.items() if diff is not None}

    for i in range(n):
        mz_i = float(features[i]["mz"])
        rt_i = float(features[i]["rt"])
        for j in range(i + 1, n):
            mz_j = float(features[j]["mz"])
            rt_j = float(features[j]["rt"])

            if abs(rt_i - rt_j) > rt_tolerance:
                continue

            mz_diff = mz_j - mz_i
            for name, expected_diff in active_diffs.items():
                if abs(abs(mz_diff) - abs(expected_diff)) <= mz_tolerance_da:
                    union(i, j)
                    if not annotations[i]:
                        annotations[i] = name.split(" vs ")[0] if mz_diff > 0 else name.split(" vs ")[1]
                    if not annotations[j]:
                        annotations[j] = name.split(" vs ")[1] if mz_diff > 0 else name.split(" vs ")[0]
                    break

    # Renumber groups
    root_map = {}
    group_counter = 0
    results = []
    for i in range(n):
        root = find_root(i)
        if root not in root_map:
            root_map[root] = group_counter
            group_counter += 1
        results.append({
            "feature_id": features[i].get("feature_id", str(i)),
            "mz": features[i]["mz"],
            "rt": features[i]["rt"],
            "group_id": root_map[root],
            "adduct_annotation": annotations[i],
        })

    return results


@click.command()
@click.option("--input", "input_file", required=True, help="Features TSV (mz, rt columns)")
@click.option("--rt-tolerance", type=float, default=5.0,
              help="RT tolerance in seconds for co-elution (default: 5)")
@click.option("--mz-tolerance", type=float, default=0.01,
              help="m/z tolerance in Da for adduct matching (default: 0.01)")
@click.option("--output", required=True, help="Output grouped TSV")
def main(input_file, rt_tolerance, mz_tolerance, output):
    features = []
    with open(input_file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            features.append(row)

    groups = find_adduct_groups(
        features, rt_tolerance=rt_tolerance, mz_tolerance_da=mz_tolerance
    )

    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["feature_id", "mz", "rt", "group_id", "adduct_annotation"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(groups)

    n_groups = len(set(g["group_id"] for g in groups))
    print(f"Grouped {len(features)} features into {n_groups} groups, written to {output}")


if __name__ == "__main__":
    main()
