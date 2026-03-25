"""
Neutral Loss Scanner
====================
Scan MS2 spectra for characteristic neutral losses from precursor ions.
Identifies fragment peaks that differ from the precursor by known neutral loss masses.

Features:
- Scan for user-defined neutral loss masses
- Configurable mass tolerance
- Works with mzML input files
- TSV output with scan details

Usage
-----
    python neutral_loss_scanner.py --input file.mzML --losses 97.977,162.053 --tolerance 0.02
    python neutral_loss_scanner.py --input file.mzML --losses 97.977 --tolerance 0.05 --output matches.tsv
"""

import csv
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def scan_neutral_losses(
    input_path: str,
    losses: list[float],
    tolerance: float = 0.02,
) -> list[dict]:
    """Scan MS2 spectra for characteristic neutral losses.

    Parameters
    ----------
    input_path : str
        Path to mzML file.
    losses : list[float]
        Neutral loss masses to search for (in Da).
    tolerance : float
        Mass tolerance in Da for matching.

    Returns
    -------
    list[dict]
        List of dicts with keys: scan_index, rt, precursor_mz, neutral_loss,
        fragment_mz, intensity, delta_da.
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    results = []
    for scan_idx in range(exp.getNrSpectra()):
        spec = exp.getSpectrum(scan_idx)
        if spec.getMSLevel() != 2:
            continue

        precursors = spec.getPrecursors()
        if not precursors:
            continue

        prec_mz = precursors[0].getMZ()
        rt = spec.getRT()
        mzs, intensities = spec.get_peaks()

        for loss in losses:
            expected_mz = prec_mz - loss
            if expected_mz <= 0:
                continue
            for i in range(len(mzs)):
                delta = abs(mzs[i] - expected_mz)
                if delta <= tolerance:
                    results.append({
                        "scan_index": scan_idx,
                        "rt": round(rt, 4),
                        "precursor_mz": round(prec_mz, 6),
                        "neutral_loss": round(loss, 6),
                        "fragment_mz": round(float(mzs[i]), 6),
                        "intensity": round(float(intensities[i]), 2),
                        "delta_da": round(delta, 6),
                    })

    return results


def create_synthetic_mzml(output_path: str, precursor_mz: float = 500.0, losses: list[float] | None = None) -> None:
    """Create a synthetic mzML file with known neutral loss peaks for testing.

    Parameters
    ----------
    output_path : str
        Path to write the synthetic mzML file.
    precursor_mz : float
        Precursor m/z value.
    losses : list[float] or None
        Neutral loss masses to embed as fragment peaks.
    """
    if losses is None:
        losses = [97.977, 162.053]

    exp = oms.MSExperiment()

    # Add MS1 spectrum
    ms1 = oms.MSSpectrum()
    ms1.setMSLevel(1)
    ms1.setRT(10.0)
    ms1.set_peaks(([precursor_mz], [10000.0]))
    exp.addSpectrum(ms1)

    # Add MS2 spectrum with neutral loss peaks
    ms2 = oms.MSSpectrum()
    ms2.setMSLevel(2)
    ms2.setRT(10.5)
    prec = oms.Precursor()
    prec.setMZ(precursor_mz)
    prec.setCharge(2)
    ms2.setPrecursors([prec])

    fragment_mzs = [precursor_mz - loss for loss in losses]
    fragment_mzs.append(150.0)  # unrelated peak
    fragment_ints = [5000.0] * len(losses) + [1000.0]
    ms2.set_peaks((sorted(fragment_mzs), [fragment_ints[i] for i in sorted(range(len(fragment_mzs)),
                  key=lambda k: fragment_mzs[k])]))
    exp.addSpectrum(ms2)

    oms.MzMLFile().store(output_path, exp)


def write_tsv(results: list[dict], output_path: str) -> None:
    """Write neutral loss scan results to TSV file.

    Parameters
    ----------
    results : list[dict]
        List of result dictionaries.
    output_path : str
        Path to output TSV file.
    """
    fieldnames = ["scan_index", "rt", "precursor_mz", "neutral_loss", "fragment_mz", "intensity", "delta_da"]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


@click.command()
@click.option("--input", "input_file", required=True, help="Path to input mzML file")
@click.option("--losses", required=True, help="Comma-separated neutral loss masses in Da")
@click.option("--tolerance", type=float, default=0.02, help="Mass tolerance in Da (default: 0.02)")
@click.option("--output", default=None, help="Output TSV file path (default: print to stdout)")
def main(input_file, losses, tolerance, output):
    losses_list = [float(x.strip()) for x in losses.split(",")]
    results = scan_neutral_losses(input_file, losses_list, tolerance)

    if output:
        write_tsv(results, output)
        print(f"Wrote {len(results)} neutral loss matches to {output}")
    else:
        print("scan_index\trt\tprecursor_mz\tneutral_loss\tfragment_mz\tintensity\tdelta_da")
        for r in results:
            print(
                f"{r['scan_index']}\t{r['rt']}\t{r['precursor_mz']}\t{r['neutral_loss']}\t"
                f"{r['fragment_mz']}\t{r['intensity']}\t{r['delta_da']}"
            )


if __name__ == "__main__":
    main()
