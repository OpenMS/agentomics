"""
Acquisition Rate Analyzer
==========================
Analyze MS1/MS2 acquisition rates over time from an mzML file.

Computes scan rates, cycle times, and duty-cycle percentages.

Usage
-----
    python acquisition_rate_analyzer.py --input run.mzML --output rates.tsv
"""

import argparse
import csv
import sys

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def analyze_acquisition_rates(exp: oms.MSExperiment) -> dict:
    """Analyze acquisition rates from an MSExperiment.

    Parameters
    ----------
    exp:
        Loaded ``pyopenms.MSExperiment``.

    Returns
    -------
    dict
        Contains per-scan records and summary statistics.
    """
    spectra = exp.getSpectra()
    if not spectra:
        return {"scans": [], "summary": {"total_scans": 0}}

    scans = []
    ms1_rts = []
    ms2_rts = []
    prev_rt = None

    for spec in spectra:
        rt = spec.getRT()
        level = spec.getMSLevel()
        delta = rt - prev_rt if prev_rt is not None else 0.0
        prev_rt = rt
        scans.append({
            "rt_sec": round(rt, 4),
            "ms_level": level,
            "delta_sec": round(delta, 4),
        })
        if level == 1:
            ms1_rts.append(rt)
        elif level == 2:
            ms2_rts.append(rt)

    rt_total = ms1_rts[-1] - ms1_rts[0] if len(ms1_rts) > 1 else 0.0

    cycle_times = []
    for i in range(1, len(ms1_rts)):
        cycle_times.append(ms1_rts[i] - ms1_rts[i - 1])

    avg_cycle = sum(cycle_times) / len(cycle_times) if cycle_times else 0.0

    ms1_rate = len(ms1_rts) / (rt_total / 60.0) if rt_total > 0 else 0.0
    ms2_rate = len(ms2_rts) / (rt_total / 60.0) if rt_total > 0 else 0.0

    ms2_per_cycle = len(ms2_rts) / len(ms1_rts) if ms1_rts else 0.0

    summary = {
        "total_scans": len(spectra),
        "ms1_count": len(ms1_rts),
        "ms2_count": len(ms2_rts),
        "rt_range_sec": round(rt_total, 2),
        "ms1_rate_per_min": round(ms1_rate, 2),
        "ms2_rate_per_min": round(ms2_rate, 2),
        "avg_cycle_time_sec": round(avg_cycle, 4),
        "avg_ms2_per_cycle": round(ms2_per_cycle, 2),
    }

    return {"scans": scans, "summary": summary}


def main():
    parser = argparse.ArgumentParser(
        description="Analyze MS1/MS2 acquisition rates over time."
    )
    parser.add_argument("--input", required=True, metavar="FILE", help="Path to mzML file")
    parser.add_argument("--output", required=True, metavar="FILE", help="Output TSV file")
    args = parser.parse_args()

    exp = oms.MSExperiment()
    oms.MzMLFile().load(args.input, exp)

    result = analyze_acquisition_rates(exp)

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["rt_sec", "ms_level", "delta_sec"], delimiter="\t")
        writer.writeheader()
        writer.writerows(result["scans"])

    s = result["summary"]
    print(f"Rates written to {args.output}")
    print(f"  Total scans      : {s['total_scans']}")
    print(f"  MS1 rate         : {s['ms1_rate_per_min']:.1f} /min")
    print(f"  MS2 rate         : {s['ms2_rate_per_min']:.1f} /min")
    print(f"  Avg cycle time   : {s['avg_cycle_time_sec']:.3f} s")


if __name__ == "__main__":
    main()
