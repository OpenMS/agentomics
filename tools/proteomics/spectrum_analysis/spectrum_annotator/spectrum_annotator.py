"""
Spectrum Annotator
==================
Annotate observed MS2 spectrum peaks with theoretical fragment ion matches.
Matches observed m/z values against theoretical b/y ions generated for a peptide.

Features:
- Match observed peaks to theoretical fragment ions
- Configurable mass tolerance
- TSV output with annotation details

Usage
-----
    python spectrum_annotator.py --mz-list "100.5,200.3,300.1" --intensities "1000,500,200" \
        --sequence PEPTIDEK --charge 2 --tolerance 0.02
    python spectrum_annotator.py --mz-list "100.5,200.3" --intensities "1000,500" \
        --sequence PEPTIDEK --charge 1 --tolerance 0.05 --output annotation.tsv
"""

import csv

import click
import pyopenms as oms


def annotate_spectrum(
    mz_values: list[float],
    intensities: list[float],
    sequence: str,
    charge: int = 1,
    tolerance: float = 0.02,
) -> list[dict]:
    """Annotate observed spectrum peaks with theoretical fragment ion matches.

    Parameters
    ----------
    mz_values : list[float]
        Observed m/z values.
    intensities : list[float]
        Observed intensity values (same length as mz_values).
    sequence : str
        Peptide amino acid sequence.
    charge : int
        Charge state for theoretical spectrum generation.
    tolerance : float
        Mass tolerance in Da for matching.

    Returns
    -------
    list[dict]
        List of dicts with keys: observed_mz, intensity, matched_ion, theoretical_mz, error_da.
    """
    aa_seq = oms.AASequence.fromString(sequence)

    # Generate theoretical spectrum
    tsg = oms.TheoreticalSpectrumGenerator()
    param = tsg.getParameters()
    param.setValue("add_b_ions", "true")
    param.setValue("add_y_ions", "true")
    param.setValue("add_a_ions", "true")
    param.setValue("add_metainfo", "true")
    tsg.setParameters(param)

    theo_spec = oms.MSSpectrum()
    tsg.getSpectrum(theo_spec, aa_seq, charge, charge)

    # Build theoretical ion list with annotations
    theo_ions = []
    theo_mzs, _ = theo_spec.get_peaks()
    annotations = theo_spec.getStringDataArrays()
    for i in range(theo_spec.size()):
        ann = ""
        if annotations and len(annotations) > 0 and i < annotations[0].size():
            raw = annotations[0][i]
            ann = raw.decode() if isinstance(raw, bytes) else str(raw)
        theo_ions.append({"mz": theo_mzs[i], "annotation": ann})

    # Match observed peaks to theoretical
    results = []
    for obs_idx in range(len(mz_values)):
        obs_mz = mz_values[obs_idx]
        obs_int = intensities[obs_idx] if obs_idx < len(intensities) else 0.0

        best_match = None
        best_error = float("inf")

        for theo in theo_ions:
            error = abs(obs_mz - theo["mz"])
            if error <= tolerance and error < best_error:
                best_error = error
                best_match = theo

        results.append({
            "observed_mz": round(obs_mz, 6),
            "intensity": round(obs_int, 2),
            "matched_ion": best_match["annotation"] if best_match else "",
            "theoretical_mz": round(best_match["mz"], 6) if best_match else "",
            "error_da": round(best_error, 6) if best_match else "",
        })

    return results


def write_tsv(results: list[dict], output_path: str) -> None:
    """Write annotation results to TSV file.

    Parameters
    ----------
    results : list[dict]
        List of annotation result dictionaries.
    output_path : str
        Path to output TSV file.
    """
    fieldnames = ["observed_mz", "intensity", "matched_ion", "theoretical_mz", "error_da"]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


@click.command(help="Annotate observed MS2 spectrum peaks with theoretical fragment ion matches.")
@click.option("--mz-list", required=True, help="Comma-separated observed m/z values")
@click.option("--intensities", required=True, help="Comma-separated observed intensities")
@click.option("--sequence", required=True, help="Peptide amino acid sequence")
@click.option("--charge", type=int, default=1, help="Charge state (default: 1)")
@click.option("--tolerance", type=float, default=0.02, help="Mass tolerance in Da (default: 0.02)")
@click.option("--output", default=None, help="Output TSV file path (default: print to stdout)")
def main(mz_list, intensities, sequence, charge, tolerance, output):
    mz_values = [float(x.strip()) for x in mz_list.split(",")]
    intensities_list = [float(x.strip()) for x in intensities.split(",")]

    results = annotate_spectrum(mz_values, intensities_list, sequence, charge, tolerance)

    if output:
        write_tsv(results, output)
        print(f"Wrote {len(results)} annotations to {output}")
    else:
        print("observed_mz\tintensity\tmatched_ion\ttheoretical_mz\terror_da")
        for r in results:
            print(f"{r['observed_mz']}\t{r['intensity']}\t{r['matched_ion']}\t{r['theoretical_mz']}\t{r['error_da']}")


if __name__ == "__main__":
    main()
