"""
Spectrum Normalizer
===================
Normalize mass spectra in mzML format using pyopenms's Normalizer.

Supports two normalization methods:
- ``to_one``: scale intensities so the maximum is 1.0
- ``to_TIC``: scale intensities so their sum equals 1.0

Usage
-----
    python spectrum_normalizer.py --input run.mzML --output normalized.mzML
    python spectrum_normalizer.py --input run.mzML --output normalized.mzML --method to_TIC
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def normalize_experiment(
    input_path: str,
    output_path: str,
    method: str = "to_one",
) -> int:
    """Normalize all spectra in an mzML file.

    Parameters
    ----------
    input_path : str
        Path to the input mzML file.
    output_path : str
        Path to write the normalized mzML file.
    method : str
        Normalization method: ``"to_one"`` (max intensity = 1.0) or
        ``"to_TIC"`` (sum of intensities = 1.0). Default ``"to_one"``.

    Returns
    -------
    int
        Number of spectra processed.
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    norm = oms.Normalizer()
    params = norm.getDefaults()
    params.setValue("method", method)
    norm.setParameters(params)
    norm.filterPeakMap(exp)

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
    "--method",
    type=click.Choice(["to_one", "to_TIC"], case_sensitive=True),
    default="to_one",
    show_default=True,
    help="Normalization method.",
)
def main(input_path: str, output_path: str, method: str) -> None:
    """Normalize spectra in an mzML file."""
    n = normalize_experiment(input_path, output_path, method)
    click.echo(f"Normalized {n} spectra ({method}) -> {output_path}")


if __name__ == "__main__":
    main()
