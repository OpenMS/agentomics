"""
Noise Filter - Threshold Mower
===============================
Remove peaks below an intensity threshold from MS spectra.

Wraps pyopenms.ThresholdMower for convenient CLI and programmatic use.

Usage
-----
    python noise_filter_threshold.py --input run.mzML --output filtered.mzML --threshold 100.0
    python noise_filter_threshold.py --input run.mzML --output filtered.mzML
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def filter_threshold(
    input_path: str,
    output_path: str,
    threshold: float = 100.0,
) -> int:
    """Filter peaks below an intensity threshold (ThresholdMower).

    All peaks with intensity strictly below ``threshold`` are removed.

    Parameters
    ----------
    input_path : str
        Path to the input mzML file.
    output_path : str
        Path to the output mzML file.
    threshold : float
        Minimum intensity threshold (default 100.0). Peaks below this
        value are removed.

    Returns
    -------
    int
        Number of spectra processed.
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    tm = oms.ThresholdMower()
    params = tm.getDefaults()
    params.setValue("threshold", float(threshold))
    tm.setParameters(params)
    tm.filterPeakMap(exp)

    oms.MzMLFile().store(output_path, exp)

    return exp.size()


@click.command(help="Filter peaks below an intensity threshold (ThresholdMower).")
@click.option("--input", "input_path", required=True, help="Input mzML file")
@click.option("--output", "output_path", required=True, help="Output mzML file")
@click.option(
    "--threshold",
    type=float,
    default=100.0,
    show_default=True,
    help="Minimum intensity threshold",
)
def main(input_path: str, output_path: str, threshold: float) -> None:
    count = filter_threshold(input_path, output_path, threshold)
    print(f"Filtered {count} spectra with threshold={threshold} -> {output_path}")


if __name__ == "__main__":
    main()
