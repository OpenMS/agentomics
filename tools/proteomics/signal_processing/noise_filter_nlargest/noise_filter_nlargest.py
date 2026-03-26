"""
Noise Filter - N Largest
=========================
Keep only the N most intense peaks in each spectrum.

Wraps pyopenms.NLargest for convenient CLI and programmatic use.

Usage
-----
    python noise_filter_nlargest.py --input run.mzML --output filtered.mzML --n 100
    python noise_filter_nlargest.py --input run.mzML --output filtered.mzML
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def filter_nlargest(
    input_path: str,
    output_path: str,
    n: int = 100,
) -> int:
    """Keep only the N most intense peaks per spectrum (NLargest).

    Parameters
    ----------
    input_path : str
        Path to the input mzML file.
    output_path : str
        Path to the output mzML file.
    n : int
        Number of most intense peaks to keep per spectrum (default 100).

    Returns
    -------
    int
        Number of spectra processed.
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    nl = oms.NLargest()
    params = nl.getDefaults()
    params.setValue("n", int(n))
    nl.setParameters(params)
    nl.filterPeakMap(exp)

    oms.MzMLFile().store(output_path, exp)

    return exp.size()


@click.command(help="Keep only the N most intense peaks per spectrum (NLargest).")
@click.option("--input", "input_path", required=True, help="Input mzML file")
@click.option("--output", "output_path", required=True, help="Output mzML file")
@click.option(
    "--n",
    type=int,
    default=100,
    show_default=True,
    help="Number of most intense peaks to keep per spectrum",
)
def main(input_path: str, output_path: str, n: int) -> None:
    count = filter_nlargest(input_path, output_path, n)
    print(f"Filtered {count} spectra, keeping top {n} peaks -> {output_path}")


if __name__ == "__main__":
    main()
