"""
Peptide Spectral Match Validator
================================
Validate PSMs by recomputing theoretical fragment ions and measuring
spectrum-to-theoretical coverage/alignment.

Usage
-----
    python peptide_spectral_match_validator.py --mzml run.mzML --peptides psms.tsv --output validation.tsv
"""

import csv
import sys
from typing import List

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit(
        "pyopenms is required. Install it with:  pip install pyopenms"
    )

PROTON = 1.007276


def generate_theoretical_spectrum(sequence: str, charge: int = 1) -> oms.MSSpectrum:
    """Generate a theoretical MS2 spectrum for a peptide.

    Parameters
    ----------
    sequence:
        Peptide amino acid sequence.
    charge:
        Maximum charge state for fragments.

    Returns
    -------
    pyopenms.MSSpectrum
        Theoretical spectrum with b and y ions.
    """
    aa_seq = oms.AASequence.fromString(sequence)
    spec = oms.MSSpectrum()
    generator = oms.TheoreticalSpectrumGenerator()

    params = generator.getParameters()
    params.setValue("add_b_ions", "true")
    params.setValue("add_y_ions", "true")
    params.setValue("add_metainfo", "true")
    generator.setParameters(params)

    generator.getSpectrum(spec, aa_seq, 1, charge)
    spec.sortByPosition()
    return spec


def validate_psm(
    experimental_spec: oms.MSSpectrum,
    sequence: str,
    tolerance: float = 0.02,
    charge: int = 1,
) -> dict:
    """Validate a single PSM by computing fragment ion coverage.

    Parameters
    ----------
    experimental_spec:
        Experimental MS2 spectrum.
    sequence:
        Peptide sequence.
    tolerance:
        Fragment mass tolerance in Da.
    charge:
        Max fragment charge.

    Returns
    -------
    dict
        Validation results including matched ions, coverage, etc.
    """
    theo_spec = generate_theoretical_spectrum(sequence, charge)

    alignment = []
    aligner = oms.SpectrumAlignment()
    params = aligner.getParameters()
    params.setValue("tolerance", float(tolerance))
    params.setValue("is_relative_tolerance", "false")
    aligner.setParameters(params)

    aligner.getSpectrumAlignment(alignment, theo_spec, experimental_spec)

    n_theoretical = theo_spec.size()
    n_matched = len(alignment)
    coverage = n_matched / n_theoretical if n_theoretical > 0 else 0.0

    aa_seq = oms.AASequence.fromString(sequence)

    return {
        "sequence": sequence,
        "theoretical_ions": n_theoretical,
        "matched_ions": n_matched,
        "coverage": round(coverage, 4),
        "peptide_mass": round(aa_seq.getMonoWeight(), 6),
    }


def load_psms(psm_path: str) -> List[dict]:
    """Load PSMs from a TSV file.

    Parameters
    ----------
    psm_path:
        Path to PSM TSV file with columns: spectrum_index, sequence, charge.

    Returns
    -------
    list
        List of PSM dicts.
    """
    psms = []
    with open(psm_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            psms.append({
                "spectrum_index": int(row["spectrum_index"]),
                "sequence": row["sequence"].strip(),
                "charge": int(row.get("charge", 1)),
            })
    return psms


def validate_psms(
    exp: oms.MSExperiment,
    psms: List[dict],
    tolerance: float = 0.02,
) -> List[dict]:
    """Validate multiple PSMs against an MSExperiment.

    Parameters
    ----------
    exp:
        Loaded MSExperiment.
    psms:
        List of PSM dicts with spectrum_index, sequence, charge.
    tolerance:
        Fragment mass tolerance in Da.

    Returns
    -------
    list
        List of validation result dicts.
    """
    spectra = exp.getSpectra()
    results = []
    for psm in psms:
        idx = psm["spectrum_index"]
        if idx < 0 or idx >= len(spectra):
            results.append({
                "spectrum_index": idx,
                "sequence": psm["sequence"],
                "theoretical_ions": 0,
                "matched_ions": 0,
                "coverage": 0.0,
                "peptide_mass": 0.0,
                "status": "invalid_index",
            })
            continue

        spec = spectra[idx]
        result = validate_psm(spec, psm["sequence"], tolerance=tolerance, charge=psm.get("charge", 1))
        result["spectrum_index"] = idx
        result["status"] = "valid" if result["coverage"] > 0 else "no_match"
        results.append(result)

    return results


def write_tsv(results: List[dict], output_path: str) -> None:
    """Write validation results to TSV.

    Parameters
    ----------
    results:
        List of result dicts.
    output_path:
        Output file path.
    """
    if not results:
        return
    fieldnames = ["spectrum_index", "sequence", "theoretical_ions", "matched_ions", "coverage",
                  "peptide_mass", "status"]
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in results:
            writer.writerow(row)


@click.command(help="Validate PSMs by recomputing fragment coverage.")
@click.option("--mzml", required=True, help="Input mzML file")
@click.option("--peptides", required=True, help="PSM TSV (spectrum_index, sequence, charge)")
@click.option("--tolerance", type=float, default=0.02, help="Fragment tolerance in Da (default: 0.02)")
@click.option("--output", required=True, help="Output TSV file path")
def main(mzml, peptides, tolerance, output):
    exp = oms.MSExperiment()
    oms.MzMLFile().load(mzml, exp)

    psms = load_psms(peptides)
    print(f"Loaded {len(psms)} PSMs")

    results = validate_psms(exp, psms, tolerance=tolerance)

    valid_count = sum(1 for r in results if r["status"] == "valid")
    print(f"Validated: {valid_count}/{len(results)} PSMs with fragment matches")

    write_tsv(results, output)
    print(f"Results written to {output}")


if __name__ == "__main__":
    main()
