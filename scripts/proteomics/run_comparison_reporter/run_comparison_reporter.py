"""
Run Comparison Reporter
========================
Compare two or more mzML files and report TIC correlation, shared
precursor m/z values, and retention-time shifts between runs.

Usage
-----
    python run_comparison_reporter.py --inputs run1.mzML run2.mzML --output comparison.json
"""

import argparse
import json
import math
import sys

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def _extract_run_info(exp: oms.MSExperiment) -> dict:
    """Extract TIC profile and precursor list from an experiment."""
    tic_profile = []
    precursor_set = set()
    rt_values = []

    for spec in exp.getSpectra():
        rt = spec.getRT()
        rt_values.append(rt)
        _, intensities = spec.get_peaks()
        tic = float(intensities.sum()) if len(intensities) > 0 else 0.0

        if spec.getMSLevel() == 1:
            tic_profile.append((rt, tic))
        elif spec.getMSLevel() == 2:
            for prec in spec.getPrecursors():
                precursor_set.add(round(prec.getMZ(), 4))

    return {
        "tic_profile": tic_profile,
        "precursors": precursor_set,
        "rt_range": (min(rt_values), max(rt_values)) if rt_values else (0.0, 0.0),
    }


def pearson_correlation(xs: list[float], ys: list[float]) -> float:
    """Compute Pearson correlation coefficient between two lists.

    Parameters
    ----------
    xs, ys:
        Equal-length numeric sequences.

    Returns
    -------
    float
        Pearson r, or 0.0 if undefined.
    """
    n = min(len(xs), len(ys))
    if n == 0:
        return 0.0
    xs, ys = xs[:n], ys[:n]
    mx = sum(xs) / n
    my = sum(ys) / n
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if sx == 0 or sy == 0:
        return 0.0
    return cov / (sx * sy)


def compare_runs(exp1: oms.MSExperiment, exp2: oms.MSExperiment) -> dict:
    """Compare two MSExperiment objects.

    Parameters
    ----------
    exp1, exp2:
        Loaded ``pyopenms.MSExperiment`` instances.

    Returns
    -------
    dict
        Comparison metrics including TIC correlation, shared precursors,
        and RT shift estimate.
    """
    info1 = _extract_run_info(exp1)
    info2 = _extract_run_info(exp2)

    tic1 = [t for _, t in info1["tic_profile"]]
    tic2 = [t for _, t in info2["tic_profile"]]
    tic_corr = pearson_correlation(tic1, tic2)

    shared = info1["precursors"] & info2["precursors"]
    only1 = info1["precursors"] - info2["precursors"]
    only2 = info2["precursors"] - info1["precursors"]

    rt_shift = info1["rt_range"][0] - info2["rt_range"][0]

    return {
        "tic_correlation": round(tic_corr, 6),
        "shared_precursors": len(shared),
        "unique_to_run1": len(only1),
        "unique_to_run2": len(only2),
        "rt_shift_sec": round(rt_shift, 2),
        "run1_ms1_tic_points": len(tic1),
        "run2_ms1_tic_points": len(tic2),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Compare mzML runs: TIC correlation, shared precursors, RT shift."
    )
    parser.add_argument(
        "--inputs", nargs=2, required=True, metavar="FILE", help="Two mzML files to compare"
    )
    parser.add_argument("--output", required=True, metavar="FILE", help="Output JSON report")
    args = parser.parse_args()

    exp1 = oms.MSExperiment()
    oms.MzMLFile().load(args.inputs[0], exp1)

    exp2 = oms.MSExperiment()
    oms.MzMLFile().load(args.inputs[1], exp2)

    result = compare_runs(exp1, exp2)

    with open(args.output, "w") as fh:
        json.dump(result, fh, indent=2)

    print(f"Comparison report written to {args.output}")
    print(f"  TIC correlation  : {result['tic_correlation']}")
    print(f"  Shared precursors: {result['shared_precursors']}")


if __name__ == "__main__":
    main()
