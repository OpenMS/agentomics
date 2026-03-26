"""
Chromatogram Peak Picker
========================
Pick peaks from chromatographic data using pyopenms PeakPickerChromatogram.

Detects peaks in chromatograms stored in mzML files and writes the picked
chromatograms to a new mzML file.

Usage
-----
    python chromatogram_peak_picker.py --input run.mzML --output picked.mzML
"""

import math
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def pick_chromatogram_peaks(input_path: str, output_path: str) -> int:
    """Pick peaks from chromatograms in an mzML file.

    Parameters
    ----------
    input_path : str
        Path to input mzML file containing chromatograms.
    output_path : str
        Path to output mzML file with picked chromatograms.

    Returns
    -------
    int
        Number of chromatograms processed.
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    chromatograms = exp.getChromatograms()
    if len(chromatograms) == 0:
        # No chromatograms, just write out the experiment as-is
        oms.MzMLFile().store(output_path, exp)
        return 0

    pp = oms.PeakPickerChromatogram()
    picked_chroms = []

    for chrom in chromatograms:
        picked_chrom = oms.MSChromatogram()
        pp.pickChromatogram(chrom, picked_chrom)
        picked_chroms.append(picked_chrom)

    # Build output experiment with picked chromatograms
    exp_out = oms.MSExperiment()
    # Copy spectra if any
    for spec in exp.getSpectra():
        exp_out.addSpectrum(spec)
    exp_out.setChromatograms(picked_chroms)

    oms.MzMLFile().store(output_path, exp_out)

    return len(picked_chroms)


def create_synthetic_chromatogram_mzml(output_path: str) -> None:
    """Create a synthetic mzML with chromatograms having Gaussian RT profiles.

    Parameters
    ----------
    output_path : str
        Path to write the synthetic mzML file.
    """
    exp = oms.MSExperiment()

    # Create two chromatograms with Gaussian peaks at different RTs
    chrom_configs = [
        {"native_id": "chrom_1", "peak_rt": 300.0, "peak_intensity": 50000.0},
        {"native_id": "chrom_2", "peak_rt": 600.0, "peak_intensity": 30000.0},
    ]

    chromatograms = []
    for config in chrom_configs:
        chrom = oms.MSChromatogram()
        chrom.setNativeID(config["native_id"])

        rt_values = []
        intensity_values = []
        sigma = 10.0  # peak width in seconds

        # Dense sampling across the RT range
        for i in range(1000):
            rt = 0.0 + i * 1.0
            intensity = 50.0  # baseline
            intensity += config["peak_intensity"] * math.exp(
                -0.5 * ((rt - config["peak_rt"]) / sigma) ** 2
            )
            rt_values.append(rt)
            intensity_values.append(intensity)

        chrom.set_peaks([rt_values, intensity_values])
        chromatograms.append(chrom)

    exp.setChromatograms(chromatograms)
    oms.MzMLFile().store(output_path, exp)


@click.command(help="Pick peaks from chromatographic data.")
@click.option("--input", "input_path", required=True, help="Input mzML file")
@click.option("--output", "output_path", required=True, help="Output mzML file with picked chromatograms")
def main(input_path: str, output_path: str) -> None:
    count = pick_chromatogram_peaks(input_path, output_path)
    click.echo(f"Picked peaks in {count} chromatograms, written to {output_path}")


if __name__ == "__main__":
    main()
