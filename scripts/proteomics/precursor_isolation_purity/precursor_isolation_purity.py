"""
Precursor Isolation Purity
===========================
Estimate precursor purity for each MS2 spectrum by examining the
surrounding MS1 spectrum within the isolation window.

Purity is defined as the fraction of total intensity in the isolation
window attributable to the target precursor ion.

Usage
-----
    python precursor_isolation_purity.py --input run.mzML --output purity.tsv
"""

import argparse
import csv
import sys

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def estimate_purity(
    ms1_spec: oms.MSSpectrum,
    precursor_mz: float,
    isolation_width: float = 2.0,
) -> float:
    """Estimate precursor purity from the MS1 spectrum.

    Parameters
    ----------
    ms1_spec:
        The MS1 spectrum closest in RT to the MS2.
    precursor_mz:
        The precursor m/z value.
    isolation_width:
        Total width of the isolation window in Da.

    Returns
    -------
    float
        Purity as a fraction (0.0 to 1.0).
    """
    mzs, intensities = ms1_spec.get_peaks()
    if len(mzs) == 0:
        return 0.0

    half_width = isolation_width / 2.0
    lo = precursor_mz - half_width
    hi = precursor_mz + half_width

    total_in_window = 0.0
    target_intensity = 0.0
    closest_dist = float("inf")

    for mz, intensity in zip(mzs, intensities):
        if lo <= mz <= hi:
            total_in_window += float(intensity)
            dist = abs(mz - precursor_mz)
            if dist < closest_dist:
                closest_dist = dist
                target_intensity = float(intensity)

    if total_in_window == 0:
        return 0.0
    return target_intensity / total_in_window


def compute_all_purities(exp: oms.MSExperiment, isolation_width: float = 2.0) -> list[dict]:
    """Compute purity for all MS2 spectra in an experiment.

    Parameters
    ----------
    exp:
        Loaded ``pyopenms.MSExperiment``.
    isolation_width:
        Isolation window width in Da.

    Returns
    -------
    list[dict]
        Each dict has: scan_index, rt, precursor_mz, purity.
    """
    spectra = exp.getSpectra()
    last_ms1 = None
    results = []

    for i, spec in enumerate(spectra):
        if spec.getMSLevel() == 1:
            last_ms1 = spec
        elif spec.getMSLevel() == 2 and last_ms1 is not None:
            for prec in spec.getPrecursors():
                prec_mz = prec.getMZ()
                iso_w = isolation_width
                if prec.getIsolationWindowLowerOffset() > 0:
                    iso_w = prec.getIsolationWindowLowerOffset() + prec.getIsolationWindowUpperOffset()
                purity = estimate_purity(last_ms1, prec_mz, iso_w)
                results.append({
                    "scan_index": i,
                    "rt": round(spec.getRT(), 4),
                    "precursor_mz": round(prec_mz, 6),
                    "purity": round(purity, 4),
                })

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Estimate precursor isolation purity from mzML."
    )
    parser.add_argument("--input", required=True, metavar="FILE", help="Path to mzML file")
    parser.add_argument("--output", required=True, metavar="FILE", help="Output purity TSV")
    parser.add_argument(
        "--isolation-width", type=float, default=2.0,
        help="Default isolation window width in Da (default: 2.0)"
    )
    args = parser.parse_args()

    exp = oms.MSExperiment()
    oms.MzMLFile().load(args.input, exp)

    purities = compute_all_purities(exp, args.isolation_width)

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["scan_index", "rt", "precursor_mz", "purity"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(purities)

    if purities:
        avg = sum(p["purity"] for p in purities) / len(purities)
        print(f"Wrote {len(purities)} purity values to {args.output}")
        print(f"  Mean purity: {avg:.4f}")
    else:
        print("No MS2 spectra with precursors found.")


if __name__ == "__main__":
    main()
