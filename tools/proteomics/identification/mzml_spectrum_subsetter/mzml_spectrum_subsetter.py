"""
mzML Spectrum Subsetter
=======================
Extract specific spectra from mzML by scan number list.

Features:
- Extract spectra by scan index (0-based)
- Output subset as new mzML file
- Preserves spectrum metadata

Usage
-----
    python mzml_spectrum_subsetter.py --input run.mzML --scans 0,1,5 --output subset.mzML
"""


import click
import pyopenms as oms


def subset_spectra(
    input_path: str,
    scan_indices: list[int],
    output_path: str,
) -> int:
    """Extract specific spectra from mzML by scan index.

    Parameters
    ----------
    input_path : str
        Path to input mzML file.
    scan_indices : list[int]
        List of 0-based scan indices to extract.
    output_path : str
        Path to output mzML file.

    Returns
    -------
    int
        Number of spectra written.
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    out_exp = oms.MSExperiment()
    n_spectra = exp.getNrSpectra()
    count = 0

    for idx in sorted(scan_indices):
        if 0 <= idx < n_spectra:
            out_exp.addSpectrum(exp.getSpectrum(idx))
            count += 1

    oms.MzMLFile().store(output_path, out_exp)
    return count


def create_synthetic_mzml(output_path: str, n_scans: int = 10) -> None:
    """Create a synthetic mzML file for testing.

    Parameters
    ----------
    output_path : str
        Path to write the synthetic mzML file.
    n_scans : int
        Number of scans to generate.
    """
    exp = oms.MSExperiment()

    for i in range(n_scans):
        spec = oms.MSSpectrum()
        spec.setMSLevel(1 if i % 2 == 0 else 2)
        spec.setRT(float(i) * 5.0)
        mzs = [100.0 + i * 10, 200.0 + i * 10, 300.0 + i * 10]
        ints = [1000.0 * (i + 1), 500.0 * (i + 1), 200.0 * (i + 1)]
        spec.set_peaks((mzs, ints))

        if spec.getMSLevel() == 2:
            prec = oms.Precursor()
            prec.setMZ(500.0)
            prec.setCharge(2)
            spec.setPrecursors([prec])

        exp.addSpectrum(spec)

    oms.MzMLFile().store(output_path, exp)


@click.command(help="Extract specific spectra from mzML by scan number list.")
@click.option("--input", "input_path", required=True, help="Path to input mzML file")
@click.option("--scans", required=True, help="Comma-separated scan indices (0-based)")
@click.option("--output", required=True, help="Path to output mzML file")
def main(input_path, scans, output):
    scan_indices = [int(x.strip()) for x in scans.split(",")]
    count = subset_spectra(input_path, scan_indices, output)
    print(f"Extracted {count} spectra to {output}")


if __name__ == "__main__":
    main()
