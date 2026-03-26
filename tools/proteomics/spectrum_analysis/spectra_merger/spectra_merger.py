"""
Spectra Merger
==============
Merge MS spectra using pyopenms SpectraMerger. Supports block-wise
merging of consecutive spectra and precursor-based merging of MS2
spectra sharing the same precursor m/z.

Features:
- Block-wise merge: combine consecutive spectra in fixed-size blocks
- Precursor merge: combine MS2 spectra with matching precursor m/z
- Configurable m/z binning width and block size

Usage
-----
    python spectra_merger.py --input run.mzML --output merged.mzML --mode block --block-size 3
    python spectra_merger.py --input run.mzML --output merged.mzML --mode precursor
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def merge_spectra(
    input_path: str,
    output_path: str,
    mode: str = "block",
    block_size: int = 3,
    mz_binning_width: float = 10.0,
    mz_binning_width_unit: str = "ppm",
) -> int:
    """Merge spectra from an mzML file.

    Parameters
    ----------
    input_path : str
        Path to the input mzML file.
    output_path : str
        Path to write the merged mzML file.
    mode : str
        Merging mode: ``"block"`` for block-wise or ``"precursor"`` for
        precursor-based merging of MS2 spectra.
    block_size : int
        Number of consecutive spectra per block (only used in block mode).
    mz_binning_width : float
        Bin width for aligning peaks across spectra (default 10.0).
    mz_binning_width_unit : str
        Unit for ``mz_binning_width``: ``"ppm"`` or ``"Da"``
        (default ``"ppm"``).

    Returns
    -------
    int
        Number of spectra in the output experiment after merging.
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    merger = oms.SpectraMerger()
    params = merger.getDefaults()
    params.setValue("mz_binning_width", float(mz_binning_width))
    params.setValue("mz_binning_width_unit", mz_binning_width_unit)

    if mode == "block":
        params.setValue("block_method:rt_block_size", int(block_size))
        merger.setParameters(params)
        merger.mergeSpectraBlockWise(exp)
    elif mode == "precursor":
        params.setValue("precursor_method:mz_tolerance", 0.1)
        params.setValue("precursor_method:rt_tolerance", 1e6)
        merger.setParameters(params)
        merger.mergeSpectraPrecursors(exp)
    else:
        raise ValueError(f"Unknown mode '{mode}'. Use 'block' or 'precursor'.")

    oms.MzMLFile().store(output_path, exp)
    return exp.getNrSpectra()


@click.command(help="Merge MS spectra from an mzML file.")
@click.option("--input", "input_path", required=True, help="Input mzML file")
@click.option("--output", "output_path", required=True, help="Output merged mzML file")
@click.option(
    "--mode",
    type=click.Choice(["block", "precursor"]),
    default="block",
    help="Merging mode (default: block)",
)
@click.option(
    "--block-size",
    type=int,
    default=3,
    help="Spectra per block in block mode (default: 3)",
)
@click.option(
    "--mz-binning-width",
    type=float,
    default=10.0,
    help="m/z binning width (default: 10.0)",
)
@click.option(
    "--mz-binning-unit",
    type=click.Choice(["ppm", "Da"]),
    default="ppm",
    help="Unit for m/z binning width (default: ppm)",
)
def main(input_path, output_path, mode, block_size, mz_binning_width, mz_binning_unit):
    n_before = oms.MSExperiment()
    oms.MzMLFile().load(input_path, n_before)
    n_before = n_before.getNrSpectra()

    n_after = merge_spectra(
        input_path,
        output_path,
        mode=mode,
        block_size=block_size,
        mz_binning_width=mz_binning_width,
        mz_binning_width_unit=mz_binning_unit,
    )
    print(f"Merged {n_before} spectra into {n_after} spectra ({mode} mode)")
    print(f"Output written to {output_path}")


if __name__ == "__main__":
    main()
