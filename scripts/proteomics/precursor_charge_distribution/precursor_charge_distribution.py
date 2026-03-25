"""
Precursor Charge Distribution
==============================
Analyze charge state distribution across MS2 spectra in mzML files.

Features:
- Count precursor charge states from MS2 spectra
- Compute percentage distribution
- TSV output

Usage
-----
    python precursor_charge_distribution.py --input run.mzML
    python precursor_charge_distribution.py --input run.mzML --output charge_dist.tsv
"""

import argparse
import csv
import sys

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def analyze_charge_distribution(input_path: str) -> list[dict]:
    """Analyze precursor charge state distribution from mzML.

    Parameters
    ----------
    input_path : str
        Path to mzML file.

    Returns
    -------
    list[dict]
        List of dicts with keys: charge, count, percentage.
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    charge_counts = {}
    total_ms2 = 0

    for i in range(exp.getNrSpectra()):
        spec = exp.getSpectrum(i)
        if spec.getMSLevel() != 2:
            continue

        total_ms2 += 1
        precursors = spec.getPrecursors()
        if precursors:
            charge = precursors[0].getCharge()
            charge_counts[charge] = charge_counts.get(charge, 0) + 1

    results = []
    for charge in sorted(charge_counts.keys()):
        count = charge_counts[charge]
        pct = (count / total_ms2 * 100) if total_ms2 > 0 else 0.0
        results.append({
            "charge": charge,
            "count": count,
            "percentage": round(pct, 2),
        })

    return results


def create_synthetic_mzml(output_path: str, charge_dist: dict[int, int] | None = None) -> None:
    """Create a synthetic mzML file with known charge distribution.

    Parameters
    ----------
    output_path : str
        Path to write the synthetic mzML file.
    charge_dist : dict or None
        Charge state to count mapping. Defaults to {2: 50, 3: 30, 4: 10, 1: 10}.
    """
    if charge_dist is None:
        charge_dist = {2: 50, 3: 30, 4: 10, 1: 10}

    exp = oms.MSExperiment()

    # MS1 survey scan
    ms1 = oms.MSSpectrum()
    ms1.setMSLevel(1)
    ms1.setRT(0.0)
    ms1.set_peaks(([500.0], [10000.0]))
    exp.addSpectrum(ms1)

    scan_idx = 1
    for charge, count in charge_dist.items():
        for _ in range(count):
            ms2 = oms.MSSpectrum()
            ms2.setMSLevel(2)
            ms2.setRT(float(scan_idx))
            prec = oms.Precursor()
            prec.setMZ(500.0)
            prec.setCharge(charge)
            ms2.setPrecursors([prec])
            ms2.set_peaks(([200.0, 300.0], [1000.0, 500.0]))
            exp.addSpectrum(ms2)
            scan_idx += 1

    oms.MzMLFile().store(output_path, exp)


def write_tsv(results: list[dict], output_path: str) -> None:
    """Write charge distribution results to TSV file.

    Parameters
    ----------
    results : list[dict]
        List of charge distribution dictionaries.
    output_path : str
        Path to output TSV file.
    """
    fieldnames = ["charge", "count", "percentage"]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


def main():
    parser = argparse.ArgumentParser(
        description="Analyze charge state distribution across MS2 spectra."
    )
    parser.add_argument("--input", required=True, help="Path to input mzML file")
    parser.add_argument("--output", default=None, help="Output TSV file path")
    args = parser.parse_args()

    results = analyze_charge_distribution(args.input)

    if args.output:
        write_tsv(results, args.output)
        print(f"Wrote charge distribution to {args.output}")
    else:
        print("charge\tcount\tpercentage")
        for r in results:
            print(f"{r['charge']}\t{r['count']}\t{r['percentage']}%")


if __name__ == "__main__":
    main()
