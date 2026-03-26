"""
Gaussian Smoother
=================
Apply Gaussian smoothing to mass spectrometry data in mzML format using
pyopenms's GaussFilter.

Gaussian smoothing reduces high-frequency noise in profile spectra by
convolving each spectrum with a Gaussian kernel of configurable width.

Usage
-----
    python gaussian_smoother.py --input run.mzML --output smoothed.mzML
    python gaussian_smoother.py --input run.mzML --output smoothed.mzML --gaussian-width 0.5
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
    gaussian_width: float = 0.2,
) -> int:
    """Apply Gaussian smoothing to all spectra in an mzML file.

    Parameters
    ----------
    input_path : str
        Path to the input mzML file.
    output_path : str
        Path to write the smoothed mzML file.
    gaussian_width : float
        Width of the Gaussian kernel in Thomson (m/z units). Default 0.2.

    Returns
    -------
    int
        Number of spectra processed.
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    gf = oms.GaussFilter()
    params = gf.getDefaults()
    params.setValue("gaussian_width", gaussian_width)
    gf.setParameters(params)
    gf.filterExperiment(exp)

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
    "--gaussian-width",
    type=float,
    default=0.2,
    show_default=True,
    help="Width of the Gaussian kernel in Thomson.",
)
def main(input_path: str, output_path: str, gaussian_width: float) -> None:
    """Apply Gaussian smoothing to spectra in an mzML file."""
    n = smooth_experiment(input_path, output_path, gaussian_width)
    click.echo(f"Smoothed {n} spectra -> {output_path}")


if __name__ == "__main__":
    main()
