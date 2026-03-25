"""
MS1 Feature Intensity Tracker
==============================
Track feature intensities across multiple mzML runs. Given a list of target
features (m/z + RT), extract the maximum intensity within tolerance from each
run's MS1 spectra.

Usage
-----
    python ms1_feature_intensity_tracker.py --inputs run1.mzML run2.mzML --features targets.tsv \\
        --ppm 10 --output tracking.tsv
"""

import csv
from typing import Dict, List

import click
import pyopenms as oms


def load_features(features_path: str) -> List[dict]:
    """Load target features from a TSV file.

    Expected columns: feature_id, mz, rt (optional).

    Parameters
    ----------
    features_path:
        Path to TSV file with target features.

    Returns
    -------
    list
        List of feature dicts.
    """
    features = []
    with open(features_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            feat = {
                "feature_id": row.get("feature_id", f"F{len(features)+1}"),
                "mz": float(row["mz"]),
            }
            if "rt" in row and row["rt"].strip():
                feat["rt"] = float(row["rt"])
            else:
                feat["rt"] = None
            features.append(feat)
    return features


def extract_intensity(
    exp: oms.MSExperiment,
    target_mz: float,
    ppm: float = 10.0,
    target_rt: float = None,
    rt_tolerance: float = 30.0,
) -> float:
    """Extract maximum intensity for a target m/z from MS1 spectra.

    Parameters
    ----------
    exp:
        Loaded MSExperiment.
    target_mz:
        Target m/z value.
    ppm:
        m/z tolerance in ppm.
    target_rt:
        Optional target RT in seconds for filtering spectra.
    rt_tolerance:
        RT tolerance in seconds if target_rt is provided.

    Returns
    -------
    float
        Maximum intensity found, or 0.0 if not detected.
    """
    import numpy as np

    mz_tol = target_mz * ppm / 1e6
    max_intensity = 0.0

    for spec in exp.getSpectra():
        if spec.getMSLevel() != 1:
            continue

        if target_rt is not None:
            rt_diff = abs(spec.getRT() - target_rt)
            if rt_diff > rt_tolerance:
                continue

        mzs, intensities = spec.get_peaks()
        if len(mzs) == 0:
            continue

        # Find peaks within tolerance
        mask = np.abs(mzs - target_mz) <= mz_tol
        if np.any(mask):
            local_max = float(np.max(intensities[mask]))
            if local_max > max_intensity:
                max_intensity = local_max

    return max_intensity


def track_features(
    run_paths: List[str],
    features: List[dict],
    ppm: float = 10.0,
    rt_tolerance: float = 30.0,
) -> List[dict]:
    """Track feature intensities across multiple runs.

    Parameters
    ----------
    run_paths:
        List of mzML file paths.
    features:
        List of target feature dicts.
    ppm:
        m/z tolerance in ppm.
    rt_tolerance:
        RT tolerance in seconds.

    Returns
    -------
    list
        List of result dicts with feature info and per-run intensities.
    """
    results = []
    for feat in features:
        row = {
            "feature_id": feat["feature_id"],
            "mz": feat["mz"],
            "rt": feat["rt"] if feat["rt"] is not None else "N/A",
        }
        for run_path in run_paths:
            exp = oms.MSExperiment()
            oms.MzMLFile().load(run_path, exp)
            intensity = extract_intensity(
                exp, feat["mz"], ppm=ppm,
                target_rt=feat["rt"], rt_tolerance=rt_tolerance,
            )
            row[run_path] = round(intensity, 4)
        results.append(row)
    return results


def track_features_from_experiments(
    experiments: Dict[str, oms.MSExperiment],
    features: List[dict],
    ppm: float = 10.0,
    rt_tolerance: float = 30.0,
) -> List[dict]:
    """Track feature intensities across pre-loaded experiments.

    Parameters
    ----------
    experiments:
        Dict mapping run name to loaded MSExperiment.
    features:
        List of target feature dicts.
    ppm:
        m/z tolerance in ppm.
    rt_tolerance:
        RT tolerance in seconds.

    Returns
    -------
    list
        List of result dicts.
    """
    results = []
    for feat in features:
        row = {
            "feature_id": feat["feature_id"],
            "mz": feat["mz"],
            "rt": feat["rt"] if feat["rt"] is not None else "N/A",
        }
        for run_name, exp in experiments.items():
            intensity = extract_intensity(
                exp, feat["mz"], ppm=ppm,
                target_rt=feat["rt"], rt_tolerance=rt_tolerance,
            )
            row[run_name] = round(intensity, 4)
        results.append(row)
    return results


def write_tsv(results: List[dict], output_path: str) -> None:
    """Write tracking results to TSV.

    Parameters
    ----------
    results:
        List of result dicts.
    output_path:
        Output file path.
    """
    if not results:
        return
    fieldnames = list(results[0].keys())
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in results:
            writer.writerow(row)


@click.command(help="Track feature intensities across multiple mzML runs.")
@click.option("--inputs", multiple=True, required=True, help="Input mzML files")
@click.option("--features", required=True, help="Target features TSV (feature_id, mz, rt)")
@click.option("--ppm", type=float, default=10.0, help="m/z tolerance in ppm (default: 10)")
@click.option("--rt-tolerance", type=float, default=30.0, help="RT tolerance in sec (default: 30)")
@click.option("--output", required=True, help="Output TSV file path")
def main(inputs, features, ppm, rt_tolerance, output):
    features_data = load_features(features)
    print(f"Loaded {len(features_data)} target features")

    results = track_features(list(inputs), features_data, ppm=ppm, rt_tolerance=rt_tolerance)

    write_tsv(results, output)
    print(f"Tracking results written to {output}")


if __name__ == "__main__":
    main()
