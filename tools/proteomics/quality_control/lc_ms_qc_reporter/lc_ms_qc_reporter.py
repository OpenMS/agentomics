"""
LC-MS QC Reporter
==================
Generate comprehensive quality-control metrics from an mzML file.

Metrics include MS1/MS2 spectrum counts, TIC stability (CV%), charge-state
distribution of precursors, and retention-time coverage.

Usage
-----
    python lc_ms_qc_reporter.py --input run.mzML --output qc_report.json
"""

import json
import math

import click
import pyopenms as oms


def compute_qc_metrics(exp: oms.MSExperiment) -> dict:
    """Compute QC metrics from an MSExperiment.

    Parameters
    ----------
    exp:
        Loaded ``pyopenms.MSExperiment`` instance.

    Returns
    -------
    dict
        Dictionary of QC metrics.
    """
    spectra = exp.getSpectra()
    ms1_count = 0
    ms2_count = 0
    tic_values = []
    charge_counts: dict[int, int] = {}
    rt_values = []

    for spec in spectra:
        level = spec.getMSLevel()
        rt = spec.getRT()
        rt_values.append(rt)
        _, intensities = spec.get_peaks()
        tic = float(intensities.sum()) if len(intensities) > 0 else 0.0

        if level == 1:
            ms1_count += 1
            tic_values.append(tic)
        elif level == 2:
            ms2_count += 1
            precursors = spec.getPrecursors()
            for prec in precursors:
                charge = prec.getCharge()
                charge_counts[charge] = charge_counts.get(charge, 0) + 1

    tic_mean = sum(tic_values) / len(tic_values) if tic_values else 0.0
    tic_std = (
        math.sqrt(sum((t - tic_mean) ** 2 for t in tic_values) / len(tic_values))
        if tic_values
        else 0.0
    )
    tic_cv = (tic_std / tic_mean * 100.0) if tic_mean > 0 else 0.0

    rt_range = (min(rt_values), max(rt_values)) if rt_values else (0.0, 0.0)

    return {
        "ms1_count": ms1_count,
        "ms2_count": ms2_count,
        "total_spectra": len(spectra),
        "tic_mean": tic_mean,
        "tic_std": tic_std,
        "tic_cv_percent": tic_cv,
        "charge_distribution": {str(k): v for k, v in sorted(charge_counts.items())},
        "rt_range_sec": list(rt_range),
    }


@click.command(help="Generate comprehensive QC report from an mzML file.")
@click.option("--input", "input", required=True, help="Path to mzML file")
@click.option("--output", required=True, help="Output JSON report path")
def main(input, output):
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input, exp)

    metrics = compute_qc_metrics(exp)

    with open(output, "w") as fh:
        json.dump(metrics, fh, indent=2)

    print(f"QC report written to {output}")
    print(f"  MS1 spectra : {metrics['ms1_count']}")
    print(f"  MS2 spectra : {metrics['ms2_count']}")
    print(f"  TIC CV%     : {metrics['tic_cv_percent']:.2f}")


if __name__ == "__main__":
    main()
