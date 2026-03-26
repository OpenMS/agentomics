"""
Savitzky-Golay Smoother
=======================
Apply Savitzky-Golay smoothing to mass spectrometry data in mzML format
using pyopenms's SavitzkyGolayFilter.

Savitzky-Golay filtering fits successive sub-sets of adjacent data points
with a low-degree polynomial by least squares, preserving peak shape better
than simple moving average filters.

Usage
-----
    python savitzky_golay_smoother.py --input run.mzML --output smoothed.mzML
    python savitzky_golay_smoother.py --input run.mzML --output smoothed.mzML --frame-length 15 --polynomial-order 4
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def smooth_experiment(
    input_path: str,
    output_path: str,
    frame_length: int = 11,
    polynomial_order: int = 3,
) -> int:
    """Apply Savitzky-Golay smoothing to all spectra in an mzML file.

    Parameters
    ----------
    input_path : str
        Path to the input mzML file.
    output_path : str
        Path to write the smoothed mzML file.
    frame_length : int
        Number of data points in the smoothing frame. Must be odd and
        greater than ``polynomial_order``. Default 11.
    polynomial_order : int
        Order of the fitting polynomial. Default 3.

    Returns
    -------
    int
        Number of spectra processed.
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    sgf = oms.SavitzkyGolayFilter()
    params = sgf.getDefaults()
    params.setValue("frame_length", frame_length)
    params.setValue("polynomial_order", polynomial_order)
    sgf.setParameters(params)
    sgf.filterExperiment(exp)

    oms.MzMLFile().store(output_path, exp)
    return exp.size()


@click.command()
@click.option(
    "--input", "input_path", required=True, help="Input mzML file path."
)
@click.option(
    "--output", "output_path", required=True, help="Output mzML file path."
)
@click.option(
    "--frame-length",
    type=int,
    default=11,
    show_default=True,
    help="Number of data points in the smoothing frame (must be odd).",
)
@click.option(
    "--polynomial-order",
    type=int,
    default=3,
    show_default=True,
    help="Order of the fitting polynomial.",
)
def main(
    input_path: str,
    output_path: str,
    frame_length: int,
    polynomial_order: int,
) -> None:
    """Apply Savitzky-Golay smoothing to spectra in an mzML file."""
    n = smooth_experiment(input_path, output_path, frame_length, polynomial_order)
    click.echo(f"Smoothed {n} spectra -> {output_path}")


if __name__ == "__main__":
    main()
