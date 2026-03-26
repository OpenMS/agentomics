"""
Signal-to-Noise Estimator
=========================
Estimate signal-to-noise ratios for peaks in MS spectra using the median method.

Uses pyopenms SignalToNoiseEstimatorMedian to compute S/N values for each peak
in every spectrum of an mzML file, and writes a TSV report.

Usage
-----
    python signal_to_noise_estimator.py --input run.mzML --output sn_report.tsv
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def estimate_sn(input_path: str, output_path: str) -> int:
    """Estimate signal-to-noise ratios for all spectra in an mzML file.

    Parameters
    ----------
    input_path : str
        Path to input mzML file.
    output_path : str
        Path to output TSV report.

    Returns
    -------
    int
        Number of spectra analyzed.
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    spectra_count = 0
    rows = []

    for spec_idx, spectrum in enumerate(exp.getSpectra()):
        if spectrum.size() == 0:
            continue

        sne = oms.SignalToNoiseEstimatorMedian()
        sne.init(spectrum)

        for i in range(spectrum.size()):
            mz = spectrum[i].getMZ()
            intensity = spectrum[i].getIntensity()
            sn_ratio = sne.getSignalToNoise(i)
            rows.append((spec_idx, mz, intensity, sn_ratio))

        spectra_count += 1

    with open(output_path, "w") as f:
        f.write("spectrum_index\tmz\tintensity\tsn_ratio\n")
        for row in rows:
            f.write(f"{row[0]}\t{row[1]:.6f}\t{row[2]:.2f}\t{row[3]:.4f}\n")

    return spectra_count


def create_synthetic_mzml_with_noise(output_path: str) -> None:
    """Create a synthetic mzML with a strong peak among noise for testing.

    Parameters
    ----------
    output_path : str
        Path to write the synthetic mzML file.
    """
    import random

    random.seed(42)

    exp = oms.MSExperiment()
    spectrum = oms.MSSpectrum()
    spectrum.setMSLevel(1)
    spectrum.setRT(60.0)

    # Generate noise peaks across a wide m/z range
    mz_values = []
    intensity_values = []

    # Low-intensity noise across the range
    for i in range(200):
        mz = 100.0 + i * 5.0
        intensity = random.uniform(10.0, 100.0)
        mz_values.append(mz)
        intensity_values.append(intensity)

    # Insert a strong peak at m/z 500
    strong_idx = None
    for i, mz in enumerate(mz_values):
        if mz >= 500.0:
            strong_idx = i
            break
    if strong_idx is not None:
        intensity_values[strong_idx] = 50000.0

    spectrum.set_peaks([mz_values, intensity_values])
    spectrum.sortByPosition()
    exp.addSpectrum(spectrum)

    oms.MzMLFile().store(output_path, exp)


@click.command(help="Estimate signal-to-noise ratios for peaks in MS spectra.")
@click.option("--input", "input_path", required=True, help="Input mzML file")
@click.option("--output", "output_path", required=True, help="Output TSV report")
def main(input_path: str, output_path: str) -> None:
    count = estimate_sn(input_path, output_path)
    click.echo(f"Analyzed {count} spectra, report written to {output_path}")


if __name__ == "__main__":
    main()
