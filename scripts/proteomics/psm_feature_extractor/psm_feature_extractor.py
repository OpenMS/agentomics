"""
PSM Feature Extractor
=====================
Extract rescoring features from PSMs by comparing experimental spectra to
theoretical spectra generated from peptide sequences.

Usage
-----
    python psm_feature_extractor.py --mzml run.mzML --peptides psms.tsv --output features.tsv
"""

import argparse
import csv
import sys
from typing import List, Optional

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

PROTON = 1.007276


def load_mzml(input_path: str) -> oms.MSExperiment:
    """Load an mzML file."""
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)
    return exp


def load_psms_tsv(psms_path: str) -> List[dict]:
    """Load PSMs from a TSV file.

    Expected columns: sequence, charge, rt, scan_index (or mz).
    """
    psms = []
    with open(psms_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            psm = {
                "sequence": row["sequence"],
                "charge": int(row["charge"]),
                "rt": float(row["rt"]),
                "scan_index": int(row.get("scan_index", -1)) if row.get("scan_index") else -1,
                "mz": float(row.get("mz", 0.0)) if row.get("mz") else 0.0,
            }
            if psm["mz"] == 0.0:
                aa_seq = oms.AASequence.fromString(psm["sequence"])
                mass = aa_seq.getMonoWeight()
                psm["mz"] = (mass + psm["charge"] * PROTON) / psm["charge"]
            psms.append(psm)
    return psms


def find_experimental_spectrum(
    exp: oms.MSExperiment,
    scan_index: int = -1,
    precursor_mz: float = 0.0,
    rt: float = 0.0,
    mz_tolerance: float = 0.02,
    rt_tolerance: float = 30.0,
) -> Optional[oms.MSSpectrum]:
    """Find an experimental MS2 spectrum by scan index or precursor m/z + RT."""
    if 0 <= scan_index < exp.size():
        s = exp[scan_index]
        if s.getMSLevel() == 2:
            return s

    # Fallback: search by precursor mz + rt
    best = None
    best_dist = float("inf")
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
        dist = abs(spectrum.getRT() - rt)
        if dist < best_dist:
            best_dist = dist
            best = spectrum

    return best


def generate_theoretical_spectrum(sequence: str, charge: int) -> oms.MSSpectrum:
    """Generate a theoretical MS2 spectrum for a peptide sequence."""
    aa_seq = oms.AASequence.fromString(sequence)
    tsg = oms.TheoreticalSpectrumGenerator()

    params = tsg.getParameters()
    params.setValue("add_b_ions", "true")
    params.setValue("add_y_ions", "true")
    params.setValue("add_metainfo", "true")
    tsg.setParameters(params)

    theo_spectrum = oms.MSSpectrum()
    tsg.getSpectrum(theo_spectrum, aa_seq, 1, charge)
    return theo_spectrum


def compute_features(
    exp_spectrum: oms.MSSpectrum,
    theo_spectrum: oms.MSSpectrum,
    sequence: str,
    charge: int,
    tolerance: float = 0.02,
) -> dict:
    """Compute rescoring features from experimental vs theoretical spectrum comparison.

    Features include:
    - matched_ions: number of matched peaks
    - total_theoretical: total theoretical peaks
    - matched_fraction: fraction of theoretical peaks matched
    - matched_intensity_fraction: fraction of experimental intensity in matched peaks
    - delta_rt: difference between experimental and expected RT (if available)
    - precursor_mass_error: mass error of precursor
    - sequence_length: length of peptide sequence
    - charge: charge state
    - num_peaks: number of peaks in experimental spectrum
    """
    # Spectrum alignment
    alignment = []
    spa = oms.SpectrumAlignment()
    params = spa.getParameters()
    params.setValue("tolerance", float(tolerance))
    params.setValue("is_relative_tolerance", "false")
    spa.setParameters(params)

    spa.getSpectrumAlignment(alignment, theo_spectrum, exp_spectrum)

    matched_ions = len(alignment)
    total_theoretical = theo_spectrum.size()

    # Compute matched intensity fraction
    exp_mz, exp_int = exp_spectrum.get_peaks()
    total_intensity = float(sum(exp_int)) if len(exp_int) > 0 else 0.0

    matched_intensity = 0.0
    matched_exp_indices = set()
    for pair in alignment:
        exp_idx = pair[1]
        if exp_idx < len(exp_int):
            matched_intensity += float(exp_int[exp_idx])
            matched_exp_indices.add(exp_idx)

    matched_intensity_fraction = matched_intensity / total_intensity if total_intensity > 0 else 0.0
    matched_fraction = matched_ions / total_theoretical if total_theoretical > 0 else 0.0

    # Compute precursor mass error
    aa_seq = oms.AASequence.fromString(sequence)
    theo_mass = aa_seq.getMonoWeight()
    theo_mz = (theo_mass + charge * PROTON) / charge

    precursors = exp_spectrum.getPrecursors()
    precursor_mass_error = 0.0
    if precursors:
        exp_mz_prec = precursors[0].getMZ()
        precursor_mass_error = (exp_mz_prec - theo_mz) * 1e6 / theo_mz if theo_mz > 0 else 0.0

    return {
        "sequence": sequence,
        "charge": charge,
        "matched_ions": matched_ions,
        "total_theoretical": total_theoretical,
        "matched_fraction": round(matched_fraction, 4),
        "matched_intensity_fraction": round(matched_intensity_fraction, 4),
        "precursor_mass_error_ppm": round(precursor_mass_error, 4),
        "sequence_length": len(sequence),
        "num_peaks": len(exp_mz),
    }


def extract_features(
    mzml_path: str,
    psms_path: str,
    output_path: str,
    tolerance: float = 0.02,
    mz_tolerance: float = 0.02,
    rt_tolerance: float = 30.0,
) -> dict:
    """Extract rescoring features for all PSMs.

    Returns statistics about the extraction.
    """
    exp = load_mzml(mzml_path)
    psms = load_psms_tsv(psms_path)

    results = []
    matched = 0

    for psm in psms:
        exp_spectrum = find_experimental_spectrum(
            exp, psm["scan_index"], psm["mz"], psm["rt"], mz_tolerance, rt_tolerance
        )
        if exp_spectrum is None:
            continue

        theo_spectrum = generate_theoretical_spectrum(psm["sequence"], psm["charge"])
        features = compute_features(exp_spectrum, theo_spectrum, psm["sequence"], psm["charge"], tolerance)
        results.append(features)
        matched += 1

    # Write output
    if results:
        with open(output_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=results[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(results)

    return {"total_psms": len(psms), "matched": matched}


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract rescoring features from PSMs.")
    parser.add_argument("--mzml", required=True, help="Input mzML file")
    parser.add_argument("--peptides", required=True, help="Input PSMs TSV file")
    parser.add_argument("--output", required=True, help="Output features TSV file")
    parser.add_argument("--tolerance", type=float, default=0.02, help="Fragment tolerance in Da (default: 0.02)")
    parser.add_argument("--mz-tolerance", type=float, default=0.02, help="Precursor m/z tolerance (default: 0.02)")
    parser.add_argument("--rt-tolerance", type=float, default=30.0, help="RT tolerance in seconds (default: 30)")
    args = parser.parse_args()

    stats = extract_features(
        args.mzml, args.peptides, args.output, args.tolerance, args.mz_tolerance, args.rt_tolerance
    )
    print(f"Extracted features for {stats['matched']} / {stats['total_psms']} PSMs to {args.output}")


if __name__ == "__main__":
    main()
