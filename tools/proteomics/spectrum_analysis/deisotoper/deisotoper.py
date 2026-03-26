"""
Deisotoper
==========
Remove isotope peaks from MS2 spectra, keeping only the monoisotopic peak
for each isotope envelope. Optionally converts peaks to singly-charged.

Wraps pyopenms.Deisotoper.deisotopeAndSingleCharge for batch processing
of all MS2 spectra in an mzML file.

Usage
-----
    python deisotoper.py --input run.mzML --output deisotoped.mzML
    python deisotoper.py --input run.mzML --output deisotoped.mzML \
        --fragment-tolerance 0.05 --min-charge 1 --max-charge 3
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def deisotope(
    input_path: str,
    output_path: str,
    fragment_tolerance: float = 0.1,
    min_charge: int = 1,
    max_charge: int = 5,
) -> int:
    """Deisotope all MS2 spectra in an mzML file.

    Parameters
    ----------
    input_path : str
        Path to input mzML file.
    output_path : str
        Path to output mzML file with deisotoped spectra.
    fragment_tolerance : float
        Fragment mass tolerance in Da for isotope peak matching.
    min_charge : int
        Minimum charge state to consider.
    max_charge : int
        Maximum charge state to consider.

    Returns
    -------
    int
        Number of spectra processed.
    """
    exp = oms.MSExperiment()
    print(f"Loading {input_path} ...")
    oms.MzMLFile().load(input_path, exp)

    count = 0
    for i in range(exp.getNrSpectra()):
        spec = exp.getSpectrum(i)
        if spec.getMSLevel() < 2:
            continue

        oms.Deisotoper.deisotopeAndSingleCharge(
            spec,
            fragment_tolerance,
            False,          # fragment_unit_ppm
            min_charge,
            max_charge,
            False,          # keep_only_deisotoped
            3,              # min_isopeaks
            10,             # max_isopeaks
            True,           # make_single_charged
            True,           # annotate_charge
            False,          # annotate_iso_peak_count
            True,           # use_decreasing_model
            2,              # start_intensity_check
            False,          # add_up_intensity
            False,          # annotate_features
        )
        exp[i] = spec
        count += 1

    oms.MzMLFile().store(output_path, exp)
    print(f"Deisotoped {count} MS2 spectra -> {output_path}")
    return count


@click.command(help="Deisotope MS2 spectra in an mzML file.")
@click.option("--input", "input_file", required=True, help="Input mzML file")
@click.option("--output", required=True, help="Output mzML file")
@click.option(
    "--fragment-tolerance",
    type=float,
    default=0.1,
    help="Fragment mass tolerance in Da (default: 0.1)",
)
@click.option(
    "--min-charge",
    type=int,
    default=1,
    help="Minimum charge state (default: 1)",
)
@click.option(
    "--max-charge",
    type=int,
    default=5,
    help="Maximum charge state (default: 5)",
)
def main(input_file, output, fragment_tolerance, min_charge, max_charge):
    count = deisotope(input_file, output, fragment_tolerance, min_charge, max_charge)
    print(f"Processed {count} spectra.")


if __name__ == "__main__":
    main()
