"""
Spectral Library Builder
========================
Build a spectral library in MSP format from mzML + peptide identification list.

Usage
-----
    python spectral_library_builder.py --input run.mzML --peptides identified.tsv --output library.msp
"""

import csv
from typing import List, Optional

import click
import pyopenms as oms

PROTON = 1.007276


def load_mzml(input_path: str) -> oms.MSExperiment:
    """Load an mzML file."""
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)
    return exp


def load_peptides_tsv(peptides_path: str) -> List[dict]:
    """Load peptide identifications from a TSV file.

    Expected columns: sequence, charge, rt, mz (optional), score (optional).
    """
    peptides = []
    with open(peptides_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pep = {
                "sequence": row["sequence"],
                "charge": int(row["charge"]),
                "rt": float(row["rt"]),
                "mz": float(row.get("mz", 0.0)) if row.get("mz") else 0.0,
                "score": float(row.get("score", 0.0)) if row.get("score") else 0.0,
            }
            if pep["mz"] == 0.0:
                aa_seq = oms.AASequence.fromString(pep["sequence"])
                mass = aa_seq.getMonoWeight()
                pep["mz"] = (mass + pep["charge"] * PROTON) / pep["charge"]
            peptides.append(pep)
    return peptides


def find_best_ms2(
    exp: oms.MSExperiment,
    precursor_mz: float,
    rt: float,
    mz_tolerance: float = 0.01,
    rt_tolerance: float = 30.0,
) -> Optional[oms.MSSpectrum]:
    """Find the best matching MS2 spectrum (highest TIC within tolerances)."""
    best = None
    best_tic = 0.0

    for spectrum in exp:
        if spectrum.getMSLevel() != 2:
            continue
        if abs(spectrum.getRT() - rt) > rt_tolerance:
            continue
        precursors = spectrum.getPrecursors()
        if not precursors:
            continue
        if abs(precursors[0].getMZ() - precursor_mz) > mz_tolerance:
            continue

        _, intensities = spectrum.get_peaks()
        tic = float(sum(intensities)) if len(intensities) > 0 else 0.0
        if tic > best_tic:
            best_tic = tic
            best = spectrum

    return best


def spectrum_to_msp(
    sequence: str,
    charge: int,
    precursor_mz: float,
    spectrum: oms.MSSpectrum,
    score: float = 0.0,
) -> str:
    """Convert a spectrum + peptide info to an MSP library entry."""
    mz_array, intensity_array = spectrum.get_peaks()
    num_peaks = len(mz_array)

    lines = []
    lines.append(f"Name: {sequence}/{charge}")
    lines.append(f"MW: {precursor_mz:.6f}")
    lines.append(f"Comment: Charge={charge} Score={score:.4f} RT={spectrum.getRT():.2f}")
    lines.append(f"Num peaks: {num_peaks}")

    for mz, intensity in zip(mz_array, intensity_array):
        lines.append(f"{mz:.6f}\t{intensity:.4f}")

    lines.append("")
    return "\n".join(lines)


def build_library(
    mzml_path: str,
    peptides_path: str,
    output_path: str,
    mz_tolerance: float = 0.01,
    rt_tolerance: float = 30.0,
) -> dict:
    """Build a spectral library from mzML + peptide identifications.

    Returns statistics about the library building.
    """
    exp = load_mzml(mzml_path)
    peptides = load_peptides_tsv(peptides_path)

    matched = 0
    with open(output_path, "w") as fh:
        for pep in peptides:
            spectrum = find_best_ms2(exp, pep["mz"], pep["rt"], mz_tolerance, rt_tolerance)
            if spectrum is not None:
                entry = spectrum_to_msp(pep["sequence"], pep["charge"], pep["mz"], spectrum, pep["score"])
                fh.write(entry + "\n")
                matched += 1

    return {
        "total_peptides": len(peptides),
        "matched_spectra": matched,
    }


@click.command(help="Build spectral library from mzML + peptide list.")
@click.option("--input", "input", required=True, help="Input mzML file")
@click.option("--peptides", required=True, help="Input peptide identifications TSV")
@click.option("--output", required=True, help="Output MSP library file")
@click.option("--mz-tolerance", type=float, default=0.01, help="m/z tolerance in Da (default: 0.01)")
@click.option("--rt-tolerance", type=float, default=30.0, help="RT tolerance in seconds (default: 30)")
def main(input, peptides, output, mz_tolerance, rt_tolerance) -> None:
    stats = build_library(input, peptides, output, mz_tolerance, rt_tolerance)
    print(f"Built library: {stats['matched_spectra']} / {stats['total_peptides']} peptides matched to spectra")


if __name__ == "__main__":
    main()
