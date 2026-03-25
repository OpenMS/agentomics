"""
Duplicate Feature Detector
============================
Detect and merge duplicate features by m/z and RT proximity.

Features within the specified m/z (ppm) and RT (seconds) tolerances
are grouped together. The feature with the highest intensity in each
group is kept as the representative.

Usage
-----
    python duplicate_feature_detector.py --input features.tsv --mz-tolerance 10 \
        --rt-tolerance 5 --output deduplicated.tsv
"""

import csv

import click
import pyopenms as oms  # noqa: F401


def detect_duplicates(
    features: list[dict],
    mz_tolerance_ppm: float = 10.0,
    rt_tolerance: float = 5.0,
) -> list[dict]:
    """Detect and group duplicate features.

    Parameters
    ----------
    features:
        List of dicts with keys: mz, rt, intensity.
    mz_tolerance_ppm:
        m/z tolerance in ppm.
    rt_tolerance:
        RT tolerance in seconds.

    Returns
    -------
    list[dict]
        Each feature augmented with ``group_id`` and ``is_duplicate``.
    """
    n = len(features)
    group_ids = list(range(n))

    def find_root(i: int) -> int:
        while group_ids[i] != i:
            group_ids[i] = group_ids[group_ids[i]]
            i = group_ids[i]
        return i

    def union(i: int, j: int):
        ri, rj = find_root(i), find_root(j)
        if ri != rj:
            group_ids[rj] = ri

    for i in range(n):
        mz_i = float(features[i]["mz"])
        rt_i = float(features[i]["rt"])
        for j in range(i + 1, n):
            mz_j = float(features[j]["mz"])
            rt_j = float(features[j]["rt"])

            mz_tol_da = mz_i * mz_tolerance_ppm / 1e6
            if abs(mz_i - mz_j) <= mz_tol_da and abs(rt_i - rt_j) <= rt_tolerance:
                union(i, j)

    # Renumber groups and find representative per group
    root_map = {}
    group_counter = 0
    groups: dict[int, list[int]] = {}

    for i in range(n):
        root = find_root(i)
        if root not in root_map:
            root_map[root] = group_counter
            group_counter += 1
        gid = root_map[root]
        groups.setdefault(gid, []).append(i)

    # Mark representatives (highest intensity per group)
    representatives = set()
    for gid, members in groups.items():
        best = max(members, key=lambda idx: float(features[idx].get("intensity", 0)))
        representatives.add(best)

    results = []
    for i in range(n):
        root = find_root(i)
        gid = root_map[root]
        feat_copy = dict(features[i])
        feat_copy["group_id"] = gid
        feat_copy["is_duplicate"] = "false" if i in representatives else "true"
        results.append(feat_copy)

    return results


def deduplicate(features: list[dict], mz_tolerance_ppm: float = 10.0, rt_tolerance: float = 5.0) -> list[dict]:
    """Return only representative (non-duplicate) features.

    Parameters
    ----------
    features:
        List of dicts with keys: mz, rt, intensity.
    mz_tolerance_ppm:
        m/z tolerance in ppm.
    rt_tolerance:
        RT tolerance in seconds.

    Returns
    -------
    list[dict]
        Only the representative features.
    """
    annotated = detect_duplicates(features, mz_tolerance_ppm, rt_tolerance)
    return [f for f in annotated if f["is_duplicate"] == "false"]


@click.command()
@click.option("--input", "input_file", required=True, help="Features TSV")
@click.option("--mz-tolerance", type=float, default=10.0,
              help="m/z tolerance in ppm (default: 10)")
@click.option("--rt-tolerance", type=float, default=5.0,
              help="RT tolerance in seconds (default: 5)")
@click.option("--output", required=True, help="Output deduplicated TSV")
def main(input_file, mz_tolerance, rt_tolerance, output):
    features = []
    with open(input_file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            features.append(row)

    deduped = deduplicate(features, mz_tolerance, rt_tolerance)

    fieldnames = list(features[0].keys()) if features else ["mz", "rt", "intensity"]
    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        # Write without the extra keys
        for d in deduped:
            row = {k: d[k] for k in fieldnames if k in d}
            writer.writerow(row)

    removed = len(features) - len(deduped)
    print(f"Deduplication: {len(features)} input, {removed} duplicates removed, {len(deduped)} kept")
    print(f"Output written to {output}")


if __name__ == "__main__":
    main()
