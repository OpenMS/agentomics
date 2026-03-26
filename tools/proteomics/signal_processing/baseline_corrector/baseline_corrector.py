"""
Baseline Corrector
==================
Remove baseline drift from mass spectrometry data in mzML format using
pyopenms's MorphologicalFilter.

The morphological filter uses top-hat filtering (a combination of erosion
and dilation) to estimate and subtract the baseline from profile spectra.

Usage
-----
    python baseline_corrector.py --input run.mzML --output corrected.mzML
    python baseline_corrector.py --input run.mzML --output corrected.mzML --struct-element-length 5.0
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def correct_baseline(
    input_path: str,
    output_path: str,
    struct_element_length: float = 3.0,
) -> int:
    """Remove baseline from all spectra in an mzML file.

    Parameters
    ----------
    input_path : str
        Path to the input mzML file.
    output_path : str
        Path to write the baseline-corrected mzML file.
    struct_element_length : float
        Length of the structuring element in Thomson (m/z units).
        Controls the scale of features considered as baseline.
        Default 3.0.

    Returns
    -------
    int
        Number of spectra processed.
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    mf = oms.MorphologicalFilter()
    params = mf.getDefaults()
    params.setValue("struc_elem_length", struct_element_length)
    mf.setParameters(params)
    mf.filterExperiment(exp)

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
    "--struct-element-length",
    type=float,
    default=3.0,
    show_default=True,
    help="Length of the structuring element in Thomson.",
)
def main(
    input_path: str, output_path: str, struct_element_length: float
) -> None:
    """Remove baseline drift from spectra in an mzML file."""
    n = correct_baseline(input_path, output_path, struct_element_length)
    click.echo(f"Baseline-corrected {n} spectra -> {output_path}")


if __name__ == "__main__":
    main()
