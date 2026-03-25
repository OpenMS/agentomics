"""
Mass Error Distribution Analyzer
=================================
Compute precursor mass-error distributions from a peptide identification
TSV and the corresponding mzML file.

For each identified peptide, the theoretical m/z is computed with pyopenms
AASequence and compared to the observed precursor m/z in the mzML.

Usage
-----
    python mass_error_distribution_analyzer.py --input peptides.tsv --mzml run.mzML --output errors.tsv
"""

import csv
import math
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

PROTON = 1.007276


def compute_mass_errors(peptide_rows: list[dict], exp: oms.MSExperiment) -> list[dict]:
    """Compute mass errors for identified peptides against spectra.

    Parameters
    ----------
    peptide_rows:
        List of dicts with keys: sequence, charge, scan_index (or precursor_mz).
    exp:
        Loaded ``pyopenms.MSExperiment``.

    Returns
    -------
    list[dict]
        Each dict has: sequence, charge, theo_mz, obs_mz, error_da, error_ppm.
    """
    spectra = exp.getSpectra()
    results = []

    for row in peptide_rows:
        sequence = row["sequence"]
        charge = int(row["charge"])

        try:
            aa = oms.AASequence.fromString(sequence)
            theo_mass = aa.getMonoWeight()
            theo_mz = (theo_mass + charge * PROTON) / charge
        except Exception:
            continue

        if "precursor_mz" in row and row["precursor_mz"]:
            obs_mz = float(row["precursor_mz"])
        elif "scan_index" in row and row["scan_index"]:
            idx = int(row["scan_index"])
            if idx >= len(spectra):
                continue
            precs = spectra[idx].getPrecursors()
            if not precs:
                continue
            obs_mz = precs[0].getMZ()
        else:
            continue

        error_da = obs_mz - theo_mz
        error_ppm = error_da / theo_mz * 1e6

        results.append({
            "sequence": sequence,
            "charge": charge,
            "theo_mz": round(theo_mz, 6),
            "obs_mz": round(obs_mz, 6),
            "error_da": round(error_da, 6),
            "error_ppm": round(error_ppm, 4),
        })

    return results


def summarize_errors(errors: list[dict]) -> dict:
    """Compute summary statistics for mass errors.

    Parameters
    ----------
    errors:
        Output of ``compute_mass_errors``.

    Returns
    -------
    dict
        Mean, std, median for both Da and ppm errors.
    """
    if not errors:
        return {"count": 0}

    ppm_vals = sorted(e["error_ppm"] for e in errors)
    da_vals = sorted(e["error_da"] for e in errors)
    n = len(ppm_vals)

    ppm_mean = sum(ppm_vals) / n
    ppm_std = math.sqrt(sum((v - ppm_mean) ** 2 for v in ppm_vals) / n)
    ppm_median = ppm_vals[n // 2]

    da_mean = sum(da_vals) / n
    da_std = math.sqrt(sum((v - da_mean) ** 2 for v in da_vals) / n)

    return {
        "count": n,
        "ppm_mean": round(ppm_mean, 4),
        "ppm_std": round(ppm_std, 4),
        "ppm_median": round(ppm_median, 4),
        "da_mean": round(da_mean, 6),
        "da_std": round(da_std, 6),
    }


@click.command(help="Compute precursor mass error distributions.")
@click.option("--input", "input", required=True, help="Peptide TSV file")
@click.option("--mzml", required=True, help="mzML file")
@click.option("--output", required=True, help="Output errors TSV")
def main(input, mzml, output):
    rows = []
    with open(input) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)

    exp = oms.MSExperiment()
    oms.MzMLFile().load(mzml, exp)

    errors = compute_mass_errors(rows, exp)

    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["sequence", "charge", "theo_mz", "obs_mz", "error_da", "error_ppm"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(errors)

    summary = summarize_errors(errors)
    print(f"Wrote {summary['count']} mass errors to {output}")
    if summary["count"] > 0:
        print(f"  Mean error: {summary['ppm_mean']:.2f} ppm (std {summary['ppm_std']:.2f})")


if __name__ == "__main__":
    main()
