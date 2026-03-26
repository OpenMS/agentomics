"""
Linear Resampler
================
Resample spectra in an mzML file to uniformly spaced m/z values using
linear interpolation.

Uses pyopenms LinearResampler to raster all spectra in an experiment to
a regular m/z grid with a configurable spacing.

Usage
-----
    python linear_resampler.py --input run.mzML --output resampled.mzML --spacing 0.01
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def resample_experiment(input_path: str, output_path: str, spacing: float = 0.01) -> int:
    """Resample spectra to uniformly spaced m/z values.

    Parameters
    ----------
    input_path : str
        Path to input mzML file.
    output_path : str
        Path to output resampled mzML file.
    spacing : float
        M/z spacing for the resampled grid.

    Returns
    -------
    int
        Number of spectra resampled.
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    lr = oms.LinearResampler()
    params = lr.getDefaults()
    params.setValue("spacing", float(spacing))
    lr.setParameters(params)

    lr.rasterExperiment(exp)

    oms.MzMLFile().store(output_path, exp)

    return exp.getNrSpectra()


def create_synthetic_irregular_mzml(output_path: str) -> None:
    """Create a synthetic mzML with irregularly spaced m/z values for testing.

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

    # Create irregularly spaced m/z points
    mz_values = sorted([500.0 + random.uniform(0, 100) for _ in range(200)])
    intensity_values = [random.uniform(100.0, 10000.0) for _ in range(200)]

    spectrum.set_peaks([mz_values, intensity_values])
    spectrum.sortByPosition()
    exp.addSpectrum(spectrum)

    oms.MzMLFile().store(output_path, exp)


@click.command(help="Resample spectra to uniformly spaced m/z values.")
@click.option("--input", "input_path", required=True, help="Input mzML file")
@click.option("--output", "output_path", required=True, help="Output resampled mzML file")
@click.option(
    "--spacing",
    default=0.01,
    type=float,
    show_default=True,
    help="M/z spacing for the resampled grid",
)
def main(input_path: str, output_path: str, spacing: float) -> None:
    count = resample_experiment(input_path, output_path, spacing)
    click.echo(f"Resampled {count} spectra with spacing {spacing}, written to {output_path}")


if __name__ == "__main__":
    main()
