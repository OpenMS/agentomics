"""
Spectrum Scoring HyperScore
============================
Score an experimental spectrum against a theoretical spectrum using a HyperScore-like
approach. Computes matched peak count, intensity dot product, and a combined score.

Features:
- Theoretical spectrum generation from peptide sequence
- Peak matching with configurable tolerance
- HyperScore-like combined scoring
- JSON output with score details

Usage
-----
    python spectrum_scoring_hyperscore.py --mz-list "100.5,200.3,300.1" --intensities "1000,500,200" \
        --sequence PEPTIDEK --charge 2
    python spectrum_scoring_hyperscore.py --mz-list "100.5,200.3" --intensities "1000,500" \
        --sequence PEPTIDEK --charge 1 --output score.json
"""

import json
import math
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def compute_hyperscore(
    mz_values: list[float],
    intensities: list[float],
    sequence: str,
    charge: int = 1,
    tolerance: float = 0.02,
) -> dict:
    """Compute a HyperScore-like score for experimental vs theoretical spectrum.

    The score is computed as: log(n_b! * n_y! * dot_product) where n_b and n_y are
    matched b-ion and y-ion counts, and dot_product is the sum of products of matched
    intensities.

    Parameters
    ----------
    mz_values : list[float]
        Experimental m/z values.
    intensities : list[float]
        Experimental intensity values.
    sequence : str
        Peptide amino acid sequence.
    charge : int
        Charge state for theoretical spectrum generation.
    tolerance : float
        Mass tolerance in Da for peak matching.

    Returns
    -------
    dict
        Dictionary with keys: hyperscore, matched_b, matched_y, matched_total,
        dot_product, sequence, charge.
    """
    aa_seq = oms.AASequence.fromString(sequence)

    # Build experimental spectrum
    exp_spec = oms.MSSpectrum()
    exp_spec.set_peaks((mz_values, intensities))
    exp_spec.sortByPosition()

    # Generate theoretical b-ions
    tsg_b = oms.TheoreticalSpectrumGenerator()
    param_b = tsg_b.getParameters()
    param_b.setValue("add_b_ions", "true")
    param_b.setValue("add_y_ions", "false")
    param_b.setValue("add_metainfo", "true")
    tsg_b.setParameters(param_b)
    theo_b = oms.MSSpectrum()
    tsg_b.getSpectrum(theo_b, aa_seq, charge, charge)

    # Generate theoretical y-ions
    tsg_y = oms.TheoreticalSpectrumGenerator()
    param_y = tsg_y.getParameters()
    param_y.setValue("add_b_ions", "false")
    param_y.setValue("add_y_ions", "true")
    param_y.setValue("add_metainfo", "true")
    tsg_y.setParameters(param_y)
    theo_y = oms.MSSpectrum()
    tsg_y.getSpectrum(theo_y, aa_seq, charge, charge)

    # Match peaks
    aligner = oms.SpectrumAlignment()
    param_align = aligner.getParameters()
    param_align.setValue("tolerance", tolerance)
    param_align.setValue("is_relative_tolerance", "false")
    aligner.setParameters(param_align)

    alignment_b = []
    if theo_b.size() > 0:
        aligner.getSpectrumAlignment(alignment_b, exp_spec, theo_b)
    matched_b = len(alignment_b)

    alignment_y = []
    if theo_y.size() > 0:
        aligner.getSpectrumAlignment(alignment_y, exp_spec, theo_y)
    matched_y = len(alignment_y)

    # Compute dot product from matched peaks
    exp_mzs, exp_ints = exp_spec.get_peaks()
    dot_product = 0.0
    for qi, _ in alignment_b:
        dot_product += float(exp_ints[qi])
    for qi, _ in alignment_y:
        dot_product += float(exp_ints[qi])

    # HyperScore = log(n_b! * n_y! * dot_product)
    if matched_b > 0 and matched_y > 0 and dot_product > 0:
        log_score = (
            math.lgamma(matched_b + 1)
            + math.lgamma(matched_y + 1)
            + math.log(dot_product)
        )
    elif (matched_b > 0 or matched_y > 0) and dot_product > 0:
        log_score = math.lgamma(max(matched_b, matched_y) + 1) + math.log(dot_product)
    else:
        log_score = 0.0

    return {
        "hyperscore": round(log_score, 6),
        "matched_b": matched_b,
        "matched_y": matched_y,
        "matched_total": matched_b + matched_y,
        "dot_product": round(dot_product, 4),
        "sequence": sequence,
        "charge": charge,
    }


@click.command(help="Score experimental spectrum against theoretical using HyperScore.")
@click.option("--mz-list", required=True, help="Comma-separated experimental m/z values")
@click.option("--intensities", required=True, help="Comma-separated experimental intensities")
@click.option("--sequence", required=True, help="Peptide amino acid sequence")
@click.option("--charge", type=int, default=1, help="Charge state (default: 1)")
@click.option("--tolerance", type=float, default=0.02, help="Mass tolerance in Da (default: 0.02)")
@click.option("--output", default=None, help="Output JSON file path (default: print to stdout)")
def main(mz_list, intensities, sequence, charge, tolerance, output):
    mz_values = [float(x.strip()) for x in mz_list.split(",")]
    intensities_list = [float(x.strip()) for x in intensities.split(",")]

    result = compute_hyperscore(mz_values, intensities_list, sequence, charge, tolerance)

    output_json = json.dumps(result, indent=2)
    if output:
        with open(output, "w") as f:
            f.write(output_json)
        print(f"Wrote score to {output}")
    else:
        print(output_json)


if __name__ == "__main__":
    main()
