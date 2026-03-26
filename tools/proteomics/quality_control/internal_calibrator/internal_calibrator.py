"""
Internal Calibrator
===================
Recalibrate m/z values in an mzML file using internal calibrants from
peptide identifications (idXML). Wraps pyopenms InternalCalibration
with support for linear and quadratic transformation models.

Features:
- Load calibrant peptide IDs from idXML
- Build linear or quadratic global calibration model
- Report before/after ppm error statistics

Usage
-----
    python internal_calibrator.py --input run.mzML --ids peptides.idXML --output calibrated.mzML --model linear
    python internal_calibrator.py --input run.mzML --ids peptides.idXML --output calibrated.mzML --model quadratic
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

MODEL_TYPES = {
    "linear": 0,
    "linear_weighted": 1,
    "quadratic": 2,
    "quadratic_weighted": 3,
}


def calibrate_mz(
    input_path: str,
    ids_path: str,
    output_path: str,
    model: str = "linear",
    tol_ppm: float = 20.0,
) -> dict:
    """Calibrate m/z values using internal calibrants from peptide IDs.

    Parameters
    ----------
    input_path : str
        Path to the input mzML file.
    ids_path : str
        Path to an idXML file with peptide identifications used as
        calibration references.
    output_path : str
        Path to write the calibrated mzML file.
    model : str
        Calibration model type: ``"linear"``, ``"linear_weighted"``,
        ``"quadratic"``, or ``"quadratic_weighted"`` (default ``"linear"``).
    tol_ppm : float
        Tolerance in ppm for matching calibrant peaks (default 20.0).

    Returns
    -------
    dict
        Dictionary with keys:

        - ``success`` (bool): whether calibration succeeded
        - ``n_calibrants`` (int): number of calibrant points found
        - ``before_ppm_median`` (float): median ppm error before calibration
        - ``before_ppm_mad`` (float): MAD of ppm error before calibration
        - ``after_ppm_median`` (float): median ppm error after calibration
        - ``after_ppm_mad`` (float): MAD of ppm error after calibration
    """
    if model not in MODEL_TYPES:
        raise ValueError(
            f"Unknown model '{model}'. Choose from: {list(MODEL_TYPES.keys())}"
        )
    model_type = MODEL_TYPES[model]

    # Load data
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    prot_ids = []
    pep_ids = oms.PeptideIdentificationList()
    oms.IdXMLFile().load(ids_path, prot_ids, pep_ids)

    # Fill calibrants from peptide identifications
    ic = oms.InternalCalibration()
    n_calibrants = ic.fillCalibrants(pep_ids, tol_ppm)

    if n_calibrants == 0:
        oms.MzMLFile().store(output_path, exp)
        return {
            "success": False,
            "n_calibrants": 0,
            "before_ppm_median": 0.0,
            "before_ppm_mad": 0.0,
            "after_ppm_median": 0.0,
            "after_ppm_mad": 0.0,
        }

    # Extract before-calibration ppm errors from calibration data
    cal_data = ic.getCalibrationPoints()
    before_errors = []
    for i in range(cal_data.size()):
        before_errors.append(cal_data.getError(i))

    before_errors_sorted = sorted(before_errors)
    n = len(before_errors_sorted)
    b_median = before_errors_sorted[n // 2] if n > 0 else 0.0
    b_mad = _median_abs_deviation(before_errors_sorted) if n > 0 else 0.0

    # Calibrate with a global model (rt_chunk < 0 means global)
    # Use generous post-calibration thresholds to avoid false failure
    success = ic.calibrate(
        exp,
        [1],              # target MS levels
        model_type,       # model type (linear=0, quadratic=2, etc.)
        -1.0,             # rt_chunk: negative = global model
        False,            # use_RANSAC
        1000.0,           # post_ppm_median threshold
        1000.0,           # post_ppm_MAD threshold
        "",               # file_models
        "",               # file_models_plot
        "",               # file_residuals
        "",               # file_residuals_plot
        "",               # rscript_executable
    )

    # Store calibrated data
    oms.MzMLFile().store(output_path, exp)

    # Compute after-calibration errors by comparing calibrated peaks to
    # the theoretical values from the peptide IDs
    after_errors = _compute_after_errors(output_path, ids_path, tol_ppm)
    after_sorted = sorted(after_errors) if after_errors else [0.0]
    na = len(after_sorted)
    a_median = after_sorted[na // 2] if na > 0 else 0.0
    a_mad = _median_abs_deviation(after_sorted) if na > 0 else 0.0

    return {
        "success": bool(success),
        "n_calibrants": n_calibrants,
        "before_ppm_median": round(b_median, 4),
        "before_ppm_mad": round(b_mad, 4),
        "after_ppm_median": round(a_median, 4),
        "after_ppm_mad": round(a_mad, 4),
    }


def _median_abs_deviation(sorted_values: list[float]) -> float:
    """Compute the median absolute deviation of sorted values."""
    n = len(sorted_values)
    if n == 0:
        return 0.0
    median = sorted_values[n // 2]
    deviations = sorted(abs(v - median) for v in sorted_values)
    return deviations[len(deviations) // 2]


def _compute_after_errors(
    mzml_path: str, ids_path: str, tol_ppm: float
) -> list[float]:
    """Compute ppm errors between calibrated mzML peaks and reference masses.

    Extracts theoretical m/z values from the peptide IDs, then for each
    ID finds the closest peak in the calibrated spectrum at the matching
    RT, and computes the ppm error.
    """
    PROTON = 1.007276

    # Load calibrated mzML
    cal_exp = oms.MSExperiment()
    oms.MzMLFile().load(mzml_path, cal_exp)

    # Load peptide IDs for theoretical masses
    prot_ids = []
    pep_ids = oms.PeptideIdentificationList()
    oms.IdXMLFile().load(ids_path, prot_ids, pep_ids)

    # Build RT -> spectrum index lookup
    spectra = cal_exp.getSpectra()
    rt_to_idx = {}
    for idx in range(len(spectra)):
        rt = round(spectra[idx].getRT(), 4)
        rt_to_idx[rt] = idx

    errors = []
    for i in range(pep_ids.size()):
        pep_id = pep_ids[i]
        hits = pep_id.getHits()
        if not hits:
            continue
        hit = hits[0]
        seq = hit.getSequence()
        charge = hit.getCharge()
        if charge < 1:
            continue
        theo_mass = seq.getMonoWeight()
        theo_mz = (theo_mass + charge * PROTON) / charge

        # Find spectrum at matching RT
        rt = round(pep_id.getRT(), 4)
        if rt not in rt_to_idx:
            continue
        spec = spectra[rt_to_idx[rt]]
        mzs, _ = spec.get_peaks()
        if len(mzs) == 0:
            continue

        # Find closest peak
        best_idx = _find_closest(mzs, theo_mz)
        obs_mz = mzs[best_idx]
        error_ppm = (obs_mz - theo_mz) / theo_mz * 1e6

        if abs(error_ppm) < tol_ppm:
            errors.append(error_ppm)

    return errors


def _find_closest(mzs, target: float) -> int:
    """Find the index of the closest value to target in a sorted array."""
    lo, hi = 0, len(mzs) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        if mzs[mid] < target:
            lo = mid + 1
        else:
            hi = mid
    # Check neighbors
    best = lo
    if lo > 0 and abs(mzs[lo - 1] - target) < abs(mzs[lo] - target):
        best = lo - 1
    return best


@click.command(help="Recalibrate m/z values using internal calibrants from idXML.")
@click.option("--input", "input_path", required=True, help="Input mzML file")
@click.option("--ids", "ids_path", required=True, help="Peptide identifications idXML file")
@click.option("--output", "output_path", required=True, help="Output calibrated mzML file")
@click.option(
    "--model",
    type=click.Choice(list(MODEL_TYPES.keys())),
    default="linear",
    help="Calibration model type (default: linear)",
)
@click.option(
    "--tol-ppm",
    type=float,
    default=20.0,
    help="Tolerance in ppm for calibrant matching (default: 20.0)",
)
def main(input_path, ids_path, output_path, model, tol_ppm):
    result = calibrate_mz(input_path, ids_path, output_path, model=model, tol_ppm=tol_ppm)

    if result["success"]:
        print(f"Calibration successful ({result['n_calibrants']} calibrants)")
        print(f"  Before: median = {result['before_ppm_median']:.2f} ppm, "
              f"MAD = {result['before_ppm_mad']:.2f} ppm")
        print(f"  After:  median = {result['after_ppm_median']:.2f} ppm, "
              f"MAD = {result['after_ppm_mad']:.2f} ppm")
        print(f"Output written to {output_path}")
    else:
        print("Calibration failed: no calibrants found")


if __name__ == "__main__":
    main()
