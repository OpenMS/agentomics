"""
MS Data ML Exporter
===================
Export MS features from mzML files as machine-learning-ready matrices (CSV).

Extracts per-spectrum features: RT, MS level, TIC, base peak m/z,
base peak intensity, number of peaks, m/z range, and intensity statistics.

Usage
-----
    python ms_data_ml_exporter.py --input run.mzML --output ml_matrix.csv
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


def extract_features(exp: oms.MSExperiment) -> List[dict]:
    """Extract ML-ready features from each spectrum in an MSExperiment.

    Parameters
    ----------
    exp:
        Loaded MSExperiment object.

    Returns
    -------
    list
        List of feature dicts, one per spectrum.
    """
    import numpy as np

    results = []
    for i, spec in enumerate(exp.getSpectra()):
        mzs, intensities = spec.get_peaks()
        n_peaks = len(mzs)
        rt = spec.getRT()
        ms_level = spec.getMSLevel()

        if n_peaks > 0:
            tic = float(np.sum(intensities))
            base_peak_idx = int(np.argmax(intensities))
            base_peak_mz = float(mzs[base_peak_idx])
            base_peak_int = float(intensities[base_peak_idx])
            mz_min = float(np.min(mzs))
            mz_max = float(np.max(mzs))
            int_mean = float(np.mean(intensities))
            int_std = float(np.std(intensities))
            int_median = float(np.median(intensities))
        else:
            tic = 0.0
            base_peak_mz = 0.0
            base_peak_int = 0.0
            mz_min = 0.0
            mz_max = 0.0
            int_mean = 0.0
            int_std = 0.0
            int_median = 0.0

        results.append({
            "spectrum_index": i,
            "rt": round(rt, 4),
            "ms_level": ms_level,
            "n_peaks": n_peaks,
            "tic": round(tic, 4),
            "base_peak_mz": round(base_peak_mz, 6),
            "base_peak_intensity": round(base_peak_int, 4),
            "mz_min": round(mz_min, 6),
            "mz_max": round(mz_max, 6),
            "intensity_mean": round(int_mean, 4),
            "intensity_std": round(int_std, 4),
            "intensity_median": round(int_median, 4),
        })

    return results


def write_csv(records: List[dict], output_path: str) -> None:
    """Write feature records to CSV.

    Parameters
    ----------
    records:
        List of feature dicts.
    output_path:
        Output CSV file path.
    """
    if not records:
        return
    fieldnames = list(records[0].keys())
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in records:
            writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(
        description="Export MS features as ML-ready matrices."
    )
    parser.add_argument("--input", required=True, help="Input mzML file")
    parser.add_argument("--output", required=True, help="Output CSV file path")
    args = parser.parse_args()

    exp = oms.MSExperiment()
    oms.MzMLFile().load(args.input, exp)

    records = extract_features(exp)
    print(f"Extracted features for {len(records)} spectra")

    write_csv(records, args.output)
    print(f"ML matrix written to {args.output}")


if __name__ == "__main__":
    main()
