"""
Spectrum Entropy Calculator
===========================
Calculate spectral entropy for MS2 spectra in mzML files.
Spectral entropy measures the information content of a spectrum.

Features:
- Normalized spectral entropy (0 to 1)
- Filter by MS level
- TSV output with scan metadata

Usage
-----
    python spectrum_entropy_calculator.py --input run.mzML --ms-level 2
    python spectrum_entropy_calculator.py --input run.mzML --ms-level 2 --output entropy.tsv
"""

import csv
import math
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def spectral_entropy(intensities: list[float]) -> float:
    """Calculate normalized spectral entropy for a list of intensities.

    Parameters
    ----------
    intensities : list[float]
        Peak intensity values.

    Returns
    -------
    float
        Normalized spectral entropy in range [0, 1].
        Returns 0 for empty spectra or spectra with one peak.
    """
    if len(intensities) <= 1:
        return 0.0

    total = sum(intensities)
    if total == 0:
        return 0.0

    # Normalize to probability distribution
    probs = [i / total for i in intensities if i > 0]

    # Shannon entropy
    entropy = -sum(p * math.log(p) for p in probs)

    # Normalize by maximum entropy (uniform distribution)
    max_entropy = math.log(len(probs))
    if max_entropy == 0:
        return 0.0

    return entropy / max_entropy


def compute_spectrum_entropies(
    input_path: str,
    ms_level: int = 2,
) -> list[dict]:
    """Compute spectral entropy for all spectra at given MS level.

    Parameters
    ----------
    input_path : str
        Path to mzML file.
    ms_level : int
        MS level to analyze (default 2).

    Returns
    -------
    list[dict]
        List of dicts with keys: scan_index, rt, n_peaks, entropy, precursor_mz.
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    results = []
    for i in range(exp.getNrSpectra()):
        spec = exp.getSpectrum(i)
        if spec.getMSLevel() != ms_level:
            continue

        rt = spec.getRT()
        _, intensities = spec.get_peaks()

        int_list = [float(x) for x in intensities]
        entropy = spectral_entropy(int_list)

        precursor_mz = 0.0
        precursors = spec.getPrecursors()
        if precursors:
            precursor_mz = precursors[0].getMZ()

        results.append({
            "scan_index": i,
            "rt": round(rt, 4),
            "n_peaks": len(int_list),
            "entropy": round(entropy, 6),
            "precursor_mz": round(precursor_mz, 6),
        })

    return results


def create_synthetic_mzml(output_path: str, n_ms2: int = 10) -> None:
    """Create a synthetic mzML file with MS2 spectra for testing.

    Parameters
    ----------
    output_path : str
        Path to write the synthetic mzML file.
    n_ms2 : int
        Number of MS2 scans to generate.
    """
    exp = oms.MSExperiment()

    # MS1
    ms1 = oms.MSSpectrum()
    ms1.setMSLevel(1)
    ms1.setRT(0.0)
    ms1.set_peaks(([500.0], [10000.0]))
    exp.addSpectrum(ms1)

    for i in range(n_ms2):
        ms2 = oms.MSSpectrum()
        ms2.setMSLevel(2)
        ms2.setRT(float(i + 1) * 2.0)

        prec = oms.Precursor()
        prec.setMZ(500.0 + i * 10)
        prec.setCharge(2)
        ms2.setPrecursors([prec])

        # Vary number of peaks and intensity distribution
        n_peaks = 5 + i * 2
        mzs = [100.0 + j * 50 for j in range(n_peaks)]
        # Varying entropy: first scans have uniform dist, later scans dominated by one peak
        ints = [1000.0 / (j + 1) ** (i * 0.3) for j in range(n_peaks)]
        ms2.set_peaks((mzs, ints))
        exp.addSpectrum(ms2)

    oms.MzMLFile().store(output_path, exp)


def write_tsv(results: list[dict], output_path: str) -> None:
    """Write entropy results to TSV file.

    Parameters
    ----------
    results : list[dict]
        List of entropy result dictionaries.
    output_path : str
        Path to output TSV file.
    """
    fieldnames = ["scan_index", "rt", "n_peaks", "entropy", "precursor_mz"]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


@click.command(help="Calculate spectral entropy for MS2 spectra in mzML.")
@click.option("--input", "input", required=True, help="Path to input mzML file")
@click.option("--ms-level", type=int, default=2, help="MS level (default: 2)")
@click.option("--output", default=None, help="Output TSV file path")
def main(input, ms_level, output):
    results = compute_spectrum_entropies(input, ms_level)

    if output:
        write_tsv(results, output)
        print(f"Wrote {len(results)} entropy values to {output}")
    else:
        print("scan_index\trt\tn_peaks\tentropy\tprecursor_mz")
        for r in results:
            print(f"{r['scan_index']}\t{r['rt']}\t{r['n_peaks']}\t{r['entropy']}\t{r['precursor_mz']}")


if __name__ == "__main__":
    main()
