"""
Precursor Recurrence Analyzer
==============================
Analyze precursor resampling in DDA runs by grouping MS2 precursors with
similar m/z and RT values.

Usage
-----
    python precursor_recurrence_analyzer.py --input run.mzML --mz-tolerance 10 --rt-tolerance 30 --output recurrence.tsv
"""

import argparse
import csv
import sys
from typing import List

try:
    import pyopenms as oms
except ImportError:
    sys.exit(
        "pyopenms is required. Install it with:  pip install pyopenms"
    )


def extract_precursors(exp: oms.MSExperiment) -> List[dict]:
    """Extract precursor information from MS2 spectra.

    Parameters
    ----------
    exp:
        Loaded MSExperiment object.

    Returns
    -------
    list
        List of dicts with spectrum_index, rt, precursor_mz, charge.
    """
    precursors = []
    for i, spec in enumerate(exp.getSpectra()):
        if spec.getMSLevel() < 2:
            continue
        rt = spec.getRT()
        for prec in spec.getPrecursors():
            precursors.append({
                "spectrum_index": i,
                "rt": rt,
                "precursor_mz": prec.getMZ(),
                "charge": prec.getCharge(),
            })
    return precursors


def find_recurrent_precursors(
    precursors: List[dict],
    mz_tolerance_ppm: float = 10.0,
    rt_tolerance_sec: float = 30.0,
) -> List[dict]:
    """Group precursors by m/z and RT proximity to identify resampling.

    Parameters
    ----------
    precursors:
        List of precursor dicts.
    mz_tolerance_ppm:
        m/z tolerance in ppm for grouping.
    rt_tolerance_sec:
        RT tolerance in seconds for grouping.

    Returns
    -------
    list
        List of group dicts with group_id, precursor m/z, RT range, count.
    """
    if not precursors:
        return []

    # Sort by m/z
    sorted_precs = sorted(precursors, key=lambda x: x["precursor_mz"])

    # Simple greedy clustering
    assigned = [False] * len(sorted_precs)
    groups = []
    group_id = 0

    for i in range(len(sorted_precs)):
        if assigned[i]:
            continue

        cluster = [sorted_precs[i]]
        assigned[i] = True
        ref_mz = sorted_precs[i]["precursor_mz"]
        ref_rt = sorted_precs[i]["rt"]

        for j in range(i + 1, len(sorted_precs)):
            if assigned[j]:
                continue
            mz_diff_ppm = abs(sorted_precs[j]["precursor_mz"] - ref_mz) / ref_mz * 1e6
            if mz_diff_ppm > mz_tolerance_ppm:
                break  # sorted by m/z, no more matches
            rt_diff = abs(sorted_precs[j]["rt"] - ref_rt)
            if rt_diff <= rt_tolerance_sec:
                cluster.append(sorted_precs[j])
                assigned[j] = True

        group_id += 1
        rts = [p["rt"] for p in cluster]
        mzs = [p["precursor_mz"] for p in cluster]
        groups.append({
            "group_id": group_id,
            "mean_mz": round(sum(mzs) / len(mzs), 6),
            "mz_range": round(max(mzs) - min(mzs), 6),
            "rt_min": round(min(rts), 2),
            "rt_max": round(max(rts), 2),
            "rt_span": round(max(rts) - min(rts), 2),
            "count": len(cluster),
            "is_recurrent": len(cluster) > 1,
        })

    return groups


def summarize_recurrence(groups: List[dict]) -> dict:
    """Summarize recurrence statistics.

    Parameters
    ----------
    groups:
        List of group dicts.

    Returns
    -------
    dict
        Summary statistics.
    """
    total_groups = len(groups)
    recurrent = [g for g in groups if g["is_recurrent"]]
    n_recurrent = len(recurrent)
    total_resampled = sum(g["count"] for g in recurrent)
    max_count = max((g["count"] for g in groups), default=0)

    return {
        "total_precursor_groups": total_groups,
        "recurrent_groups": n_recurrent,
        "unique_groups": total_groups - n_recurrent,
        "total_resampled_scans": total_resampled,
        "max_resampling": max_count,
        "recurrence_rate": round(n_recurrent / total_groups, 4) if total_groups > 0 else 0.0,
    }


def write_tsv(groups: List[dict], output_path: str) -> None:
    """Write group results to TSV.

    Parameters
    ----------
    groups:
        List of group dicts.
    output_path:
        Output file path.
    """
    if not groups:
        return
    fieldnames = list(groups[0].keys())
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in groups:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(
        description="Analyze precursor resampling in DDA runs."
    )
    parser.add_argument("--input", required=True, help="Input mzML file")
    parser.add_argument("--mz-tolerance", type=float, default=10.0, help="m/z tolerance in ppm (default: 10)")
    parser.add_argument("--rt-tolerance", type=float, default=30.0, help="RT tolerance in seconds (default: 30)")
    parser.add_argument("--output", default=None, help="Output TSV file path")
    args = parser.parse_args()

    exp = oms.MSExperiment()
    oms.MzMLFile().load(args.input, exp)

    precursors = extract_precursors(exp)
    print(f"Extracted {len(precursors)} precursors")

    groups = find_recurrent_precursors(precursors, args.mz_tolerance, args.rt_tolerance)
    summary = summarize_recurrence(groups)

    print(f"Precursor groups  : {summary['total_precursor_groups']}")
    print(f"Recurrent groups  : {summary['recurrent_groups']}")
    print(f"Unique groups     : {summary['unique_groups']}")
    print(f"Recurrence rate   : {summary['recurrence_rate']:.2%}")
    print(f"Max resampling    : {summary['max_resampling']}")

    if args.output:
        write_tsv(groups, args.output)
        print(f"\nResults written to {args.output}")


if __name__ == "__main__":
    main()
