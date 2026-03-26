"""
Metabolite Spectral Matcher
============================
Match experimental MS/MS spectra against a spectral library for metabolite
identification.

Wraps pyopenms.MetaboliteSpectralMatching to compare query spectra against
library spectra and return match scores.

Usage
-----
    python metabolite_spectral_matcher.py --input spectra.mzML \
        --library lib.mzML --output matches.mzTab --precursor-tol 0.1
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def match_spectra(
    input_path: str,
    library_path: str,
    output_path: str,
    precursor_tol: float = 0.1,
) -> int:
    """Match query spectra against a spectral library.

    Parameters
    ----------
    input_path : str
        Path to input mzML file with query MS/MS spectra.
    library_path : str
        Path to spectral library mzML file.
    output_path : str
        Path to output mzTab file.
    precursor_tol : float
        Precursor mass tolerance in Da (default 0.1).

    Returns
    -------
    int
        Number of spectra matched.
    """
    # Load query spectra
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    # Load library spectra
    lib = oms.MSExperiment()
    oms.MzMLFile().load(library_path, lib)

    # Configure and run spectral matching
    msm = oms.MetaboliteSpectralMatching()
    params = msm.getDefaults()
    params.setValue("prec_mass_error_value", precursor_tol)
    msm.setParameters(params)

    mztab = oms.MzTab()
    out_spectra = oms.String("")
    msm.run(exp, lib, mztab, out_spectra)

    # Store results
    oms.MzTabFile().store(output_path, mztab)

    # Count matches by reading SML lines in mzTab
    match_count = 0
    with open(output_path) as f:
        for line in f:
            if line.startswith("SML\t"):
                match_count += 1

    return match_count


def create_synthetic_query_mzml(output_path: str, precursor_mz: float = 181.07) -> None:
    """Create a synthetic mzML with MS2 spectra for testing.

    Creates 3 MS2 spectra with different precursor m/z values.
    """
    exp = oms.MSExperiment()
    mzs_list = [
        [50.0, 73.0, 89.0, 119.0, 163.0],
        [55.0, 80.0, 95.0, 125.0, 170.0],
        [60.0, 85.0, 100.0, 130.0, 175.0],
    ]
    intensities = [100.0, 200.0, 500.0, 300.0, 150.0]

    for i, mzs in enumerate(mzs_list):
        spec = oms.MSSpectrum()
        spec.setMSLevel(2)
        spec.setRT(100.0 + i * 10.0)

        prec = oms.Precursor()
        prec.setMZ(precursor_mz + i * 100.0)
        prec.setCharge(1)
        spec.setPrecursors([prec])

        spec.set_peaks((mzs, intensities))
        exp.addSpectrum(spec)

    oms.MzMLFile().store(output_path, exp)


def create_synthetic_library_mzml(output_path: str, precursor_mz: float = 181.07) -> None:
    """Create a synthetic spectral library mzML for testing.

    Creates matching library spectra with the same fragment peaks.
    """
    lib = oms.MSExperiment()
    mzs_list = [
        [50.0, 73.0, 89.0, 119.0, 163.0],
        [55.0, 80.0, 95.0, 125.0, 170.0],
        [60.0, 85.0, 100.0, 130.0, 175.0],
    ]
    intensities = [100.0, 200.0, 500.0, 300.0, 150.0]

    for i, mzs in enumerate(mzs_list):
        spec = oms.MSSpectrum()
        spec.setMSLevel(2)
        spec.setRT(float(i))
        spec.setName(f"Compound_{i}")

        prec = oms.Precursor()
        prec.setMZ(precursor_mz + i * 100.0)
        prec.setCharge(1)
        spec.setPrecursors([prec])

        spec.set_peaks((mzs, intensities))
        lib.addSpectrum(spec)

    oms.MzMLFile().store(output_path, lib)


@click.command(help="Match MS/MS spectra against a spectral library.")
@click.option("--input", "input_path", required=True, help="Input mzML file with query spectra")
@click.option("--library", "library_path", required=True, help="Spectral library mzML file")
@click.option("--output", "output_path", required=True, help="Output mzTab file")
@click.option("--precursor-tol", default=0.1, type=float, help="Precursor mass tolerance in Da")
def main(input_path, library_path, output_path, precursor_tol) -> None:
    matches = match_spectra(input_path, library_path, output_path, precursor_tol)
    click.echo(f"Matched {matches} spectra, saved to {output_path}")


if __name__ == "__main__":
    main()
