"""
Noise Filter - Window Mower
============================
Remove low-intensity peaks from MS spectra using a sliding window approach.
Within each window of a given m/z width, only the N most intense peaks are kept.

Wraps pyopenms.WindowMower for convenient CLI and programmatic use.

Usage
-----
    python noise_filter_window.py --input run.mzML --output filtered.mzML --window-size 50.0 --peak-count 3
    python noise_filter_window.py --input run.mzML --output filtered.mzML
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def filter_window(
    input_path: str,
    output_path: str,
    window_size: float = 50.0,
    peak_count: int = 3,
) -> int:
    """Filter peaks using a sliding window approach (WindowMower).

    Within each window of ``window_size`` Da, only the ``peak_count`` most
    intense peaks are retained.

    Parameters
    ----------
    input_path : str
        Path to the input mzML file.
    output_path : str
        Path to the output mzML file.
    window_size : float
        Width of the sliding window in Da (default 50.0).
    peak_count : int
        Number of peaks to keep per window (default 3).

    Returns
    -------
    int
        Number of spectra processed.
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    wm = oms.WindowMower()
    params = wm.getDefaults()
    params.setValue("windowsize", float(window_size))
    params.setValue("peakcount", int(peak_count))
    wm.setParameters(params)
    wm.filterPeakMap(exp)

    oms.MzMLFile().store(output_path, exp)

    return exp.size()


@click.command(help="Filter peaks using a sliding window approach (WindowMower).")
@click.option("--input", "input_path", required=True, help="Input mzML file")
@click.option("--output", "output_path", required=True, help="Output mzML file")
@click.option(
    "--window-size",
    type=float,
    default=50.0,
    show_default=True,
    help="Window size in Da",
)
@click.option(
    "--peak-count",
    type=int,
    default=3,
    show_default=True,
    help="Number of peaks to keep per window",
)
def main(input_path: str, output_path: str, window_size: float, peak_count: int) -> None:
    count = filter_window(input_path, output_path, window_size, peak_count)
    print(
        f"Filtered {count} spectra with window_size={window_size}, "
        f"peak_count={peak_count} -> {output_path}"
    )


if __name__ == "__main__":
    main()
