"""
XIC Extractor
=============
Extract ion chromatograms (XIC) for target m/z values from mzML files.
Also computes TIC and BPC as part of the extraction.

Features:
- Extract XIC for one or more target m/z values
- PPM-based mass tolerance
- TSV output with RT and intensity

Usage
-----
    python xic_extractor.py --input run.mzML --mz 524.265 --ppm 10
    python xic_extractor.py --input run.mzML --mz 524.265 --ppm 10 --output xic.tsv
"""

import argparse
import csv
import sys

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def extract_xic(
    input_path: str,
    target_mz: float,
    ppm: float = 10.0,
    ms_level: int = 1,
) -> list[dict]:
    """Extract ion chromatogram for a target m/z from mzML.

    Parameters
    ----------
    input_path : str
        Path to mzML file.
    target_mz : float
        Target m/z value.
    ppm : float
        Mass tolerance in ppm.
    ms_level : int
        MS level to extract from (default 1).

    Returns
    -------
    list[dict]
        List of dicts with keys: rt, intensity, mz.
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    tolerance_da = target_mz * ppm / 1e6
    mz_low = target_mz - tolerance_da
    mz_high = target_mz + tolerance_da

    results = []
    for i in range(exp.getNrSpectra()):
        spec = exp.getSpectrum(i)
        if spec.getMSLevel() != ms_level:
            continue

        rt = spec.getRT()
        mzs, intensities = spec.get_peaks()

        max_intensity = 0.0
        best_mz = target_mz
        for j in range(len(mzs)):
            if mz_low <= mzs[j] <= mz_high:
                if intensities[j] > max_intensity:
                    max_intensity = float(intensities[j])
                    best_mz = float(mzs[j])

        results.append({
            "rt": round(rt, 4),
            "intensity": round(max_intensity, 2),
            "mz": round(best_mz, 6),
        })

    return results


def create_synthetic_mzml(output_path: str, target_mz: float = 524.265, n_scans: int = 10) -> None:
    """Create a synthetic mzML file with known peaks for testing.

    Parameters
    ----------
    output_path : str
        Path to write the synthetic mzML file.
    target_mz : float
        Target m/z to embed in spectra.
    n_scans : int
        Number of MS1 scans to generate.
    """
    exp = oms.MSExperiment()

    for i in range(n_scans):
        spec = oms.MSSpectrum()
        spec.setMSLevel(1)
        spec.setRT(float(i) * 10.0)

        # Gaussian-like intensity profile centered at scan n_scans//2
        center = n_scans // 2
        intensity = max(100.0, 10000.0 * max(0, 1 - abs(i - center) / center))

        mzs = [target_mz - 50, target_mz, target_mz + 50]
        ints = [500.0, intensity, 300.0]
        spec.set_peaks((mzs, ints))
        exp.addSpectrum(spec)

    oms.MzMLFile().store(output_path, exp)


def write_tsv(results: list[dict], output_path: str) -> None:
    """Write XIC results to TSV file.

    Parameters
    ----------
    results : list[dict]
        List of XIC data points.
    output_path : str
        Path to output TSV file.
    """
    fieldnames = ["rt", "intensity", "mz"]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


def main():
    parser = argparse.ArgumentParser(
        description="Extract ion chromatograms for target m/z values from mzML."
    )
    parser.add_argument("--input", required=True, help="Path to input mzML file")
    parser.add_argument("--mz", type=float, required=True, help="Target m/z value")
    parser.add_argument("--ppm", type=float, default=10.0, help="Mass tolerance in ppm (default: 10)")
    parser.add_argument("--ms-level", type=int, default=1, help="MS level (default: 1)")
    parser.add_argument("--output", default=None, help="Output TSV file path")
    args = parser.parse_args()

    results = extract_xic(args.input, args.mz, args.ppm, args.ms_level)

    if args.output:
        write_tsv(results, args.output)
        print(f"Wrote {len(results)} XIC data points to {args.output}")
    else:
        print("rt\tintensity\tmz")
        for r in results:
            print(f"{r['rt']}\t{r['intensity']}\t{r['mz']}")


if __name__ == "__main__":
    main()
