"""
Peak Picker HiRes
=================
Pick peaks from profile (raw) mass spectra using the high-resolution peak picker.

Converts profile-mode mzML data to centroided mzML by detecting peaks in the
profile data using pyopenms PeakPickerHiRes.

Usage
-----
    python peak_picker_hires.py --input profile.mzML --output centroid.mzML --signal-to-noise 1.0
"""

import math
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def pick_peaks(input_path: str, output_path: str, signal_to_noise: float = 1.0) -> int:
    """Pick peaks from profile spectra using PeakPickerHiRes.

    Parameters
    ----------
    input_path : str
        Path to input profile mzML file.
    output_path : str
        Path to output centroided mzML file.
    signal_to_noise : float
        Minimum signal-to-noise ratio for a peak to be picked.

    Returns
    -------
    int
        Number of spectra processed.
    """
    exp_in = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp_in)

    exp_out = oms.MSExperiment()

    pp = oms.PeakPickerHiRes()
    params = pp.getDefaults()
    params.setValue("signal_to_noise", float(signal_to_noise))
    pp.setParameters(params)

    for spectrum in exp_in.getSpectra():
        picked = oms.MSSpectrum()
        pp.pick(spectrum, picked)
        exp_out.addSpectrum(picked)

    oms.MzMLFile().store(output_path, exp_out)

    return exp_out.getNrSpectra()


def create_synthetic_profile_mzml(output_path: str) -> None:
    """Create a synthetic profile-mode mzML with Gaussian-shaped peaks for testing.

    Parameters
    ----------
    output_path : str
        Path to write the synthetic mzML file.
    """
    exp = oms.MSExperiment()
    spectrum = oms.MSSpectrum()
    spectrum.setMSLevel(1)
    spectrum.setRT(60.0)

    mz_values = []
    intensity_values = []

    # Create dense m/z points with Gaussian peaks at 400, 500, and 600 m/z.
    # Spacing must be much smaller than sigma for the peak picker to work.
    peak_centers = [400.0, 500.0, 600.0]
    peak_intensities = [10000.0, 50000.0, 25000.0]
    sigma = 0.05  # peak width in Da

    # Dense sampling from 350 to 650: 0.01 Da spacing gives ~10 points per sigma
    n_points = 30000
    for i in range(n_points):
        mz = 350.0 + i * (300.0 / n_points)
        intensity = 10.0  # baseline noise
        for center, height in zip(peak_centers, peak_intensities):
            intensity += height * math.exp(-0.5 * ((mz - center) / sigma) ** 2)
        mz_values.append(mz)
        intensity_values.append(intensity)

    spectrum.set_peaks([mz_values, intensity_values])
    spectrum.sortByPosition()
    exp.addSpectrum(spectrum)

    oms.MzMLFile().store(output_path, exp)


@click.command(help="Pick peaks from profile mass spectra.")
@click.option("--input", "input_path", required=True, help="Input profile mzML file")
@click.option("--output", "output_path", required=True, help="Output centroided mzML file")
@click.option(
    "--signal-to-noise",
    default=1.0,
    type=float,
    show_default=True,
    help="Minimum signal-to-noise ratio",
)
def main(input_path: str, output_path: str, signal_to_noise: float) -> None:
    count = pick_peaks(input_path, output_path, signal_to_noise)
    click.echo(f"Picked peaks in {count} spectra, written to {output_path}")


if __name__ == "__main__":
    main()
