"""
Sample Complexity Estimator
============================
Estimate sample complexity from MS1 peak density in mzML files.

Reports peak density per RT window, total unique peak count, and
a complexity score based on average peak density.

Usage
-----
    python sample_complexity_estimator.py --input run.mzML --output complexity.json
"""

import json
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit(
        "pyopenms is required. Install it with:  pip install pyopenms"
    )


def estimate_complexity(
    exp: oms.MSExperiment,
    intensity_threshold: float = 0.0,
) -> dict:
    """Estimate sample complexity from MS1 peak density.

    Parameters
    ----------
    exp:
        Loaded MSExperiment object.
    intensity_threshold:
        Minimum intensity to count a peak.

    Returns
    -------
    dict
        Complexity metrics.
    """
    import numpy as np

    ms1_spectra = [s for s in exp.getSpectra() if s.getMSLevel() == 1]
    if not ms1_spectra:
        return {
            "n_ms1_spectra": 0,
            "total_peaks": 0,
            "avg_peaks_per_spectrum": 0.0,
            "max_peaks_per_spectrum": 0,
            "min_peaks_per_spectrum": 0,
            "complexity_score": "N/A",
            "per_spectrum": [],
        }

    per_spectrum = []
    peak_counts = []

    for spec in ms1_spectra:
        mzs, intensities = spec.get_peaks()
        if intensity_threshold > 0 and len(intensities) > 0:
            mask = intensities >= intensity_threshold
            n_peaks = int(np.sum(mask))
        else:
            n_peaks = len(mzs)

        rt = spec.getRT()
        per_spectrum.append({
            "rt": round(rt, 4),
            "n_peaks": n_peaks,
        })
        peak_counts.append(n_peaks)

    total_peaks = sum(peak_counts)
    avg_peaks = total_peaks / len(peak_counts)
    max_peaks = max(peak_counts)
    min_peaks = min(peak_counts)

    # Simple complexity classification
    if avg_peaks > 5000:
        complexity = "very_high"
    elif avg_peaks > 1000:
        complexity = "high"
    elif avg_peaks > 200:
        complexity = "medium"
    elif avg_peaks > 50:
        complexity = "low"
    else:
        complexity = "very_low"

    return {
        "n_ms1_spectra": len(ms1_spectra),
        "total_peaks": total_peaks,
        "avg_peaks_per_spectrum": round(avg_peaks, 2),
        "max_peaks_per_spectrum": max_peaks,
        "min_peaks_per_spectrum": min_peaks,
        "complexity_score": complexity,
        "per_spectrum": per_spectrum,
    }


def write_json(result: dict, output_path: str) -> None:
    """Write complexity result to JSON.

    Parameters
    ----------
    result:
        Complexity metrics dict.
    output_path:
        Output file path.
    """
    with open(output_path, "w") as fh:
        json.dump(result, fh, indent=2)


@click.command(help="Estimate sample complexity from MS1 peak density.")
@click.option("--input", "input", required=True, help="Input mzML file")
@click.option("--intensity-threshold", type=float, default=0.0, help="Minimum intensity to count a peak (default: 0)")
@click.option("--output", default=None, help="Output JSON file path")
def main(input, intensity_threshold, output):
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input, exp)

    result = estimate_complexity(exp, intensity_threshold=intensity_threshold)

    print(f"MS1 spectra       : {result['n_ms1_spectra']}")
    print(f"Total peaks       : {result['total_peaks']}")
    print(f"Avg peaks/spectrum: {result['avg_peaks_per_spectrum']}")
    print(f"Max peaks/spectrum: {result['max_peaks_per_spectrum']}")
    print(f"Complexity score  : {result['complexity_score']}")

    if output:
        write_json(result, output)
        print(f"\nResults written to {output}")


if __name__ == "__main__":
    main()
