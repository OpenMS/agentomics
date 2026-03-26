"""
Peak Type Estimator
===================
Estimate whether spectra in an mzML file contain profile or centroided data.

Uses pyopenms PeakTypeEstimator to classify each spectrum and produces a
TSV report with the results.

Usage
-----
    python peak_type_estimator.py --input run.mzML --output report.tsv
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

# Map pyopenms SpectrumSettings.SpectrumType enum values to human-readable names
_PEAK_TYPE_NAMES = {
    oms.SpectrumSettings.SpectrumType.UNKNOWN: "unknown",
    oms.SpectrumSettings.SpectrumType.CENTROID: "centroid",
    oms.SpectrumSettings.SpectrumType.PROFILE: "profile",
}


def estimate_peak_type(input_path: str, output_path: str) -> dict:
    """Estimate peak type (profile/centroid) for each spectrum in an mzML file.

    Parameters
    ----------
    input_path : str
        Path to input mzML file.
    output_path : str
        Path to output TSV report.

    Returns
    -------
    dict
        Dictionary with keys: 'type_counts' (mapping type name to count),
        'total_spectra' (int).
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    pte = oms.PeakTypeEstimator()

    type_counts = {}
    rows = []

    for spec_idx, spectrum in enumerate(exp.getSpectra()):
        ms_level = spectrum.getMSLevel()
        peak_type_enum = pte.estimateType(spectrum)
        peak_type_name = _PEAK_TYPE_NAMES.get(peak_type_enum, "unknown")

        rows.append((spec_idx, ms_level, peak_type_name))
        type_counts[peak_type_name] = type_counts.get(peak_type_name, 0) + 1

    with open(output_path, "w") as f:
        f.write("spectrum_index\tms_level\tpeak_type\n")
        for row in rows:
            f.write(f"{row[0]}\t{row[1]}\t{row[2]}\n")

    return {
        "type_counts": type_counts,
        "total_spectra": len(rows),
    }


def create_synthetic_profile_spectrum() -> oms.MSSpectrum:
    """Create a synthetic profile spectrum with dense m/z points.

    Returns
    -------
    oms.MSSpectrum
        A profile-mode spectrum.
    """
    import math

    spectrum = oms.MSSpectrum()
    spectrum.setMSLevel(1)
    spectrum.setRT(60.0)

    mz_values = []
    intensity_values = []

    # Dense sampling with broad Gaussian peaks spanning the full range.
    # PeakTypeEstimator needs many points with broad intensity variation
    # to reliably classify as profile.
    n_points = 10000
    for i in range(n_points):
        mz = 400.0 + i * 0.01
        intensity = 50.0
        # Broad overlapping peaks across the full m/z range
        for center in range(410, 500, 5):
            intensity += 3000.0 * math.exp(-0.5 * ((mz - center) / 1.0) ** 2)
        mz_values.append(mz)
        intensity_values.append(intensity)

    spectrum.set_peaks([mz_values, intensity_values])
    return spectrum


def create_synthetic_centroid_spectrum() -> oms.MSSpectrum:
    """Create a synthetic centroided spectrum with sparse m/z points.

    Returns
    -------
    oms.MSSpectrum
        A centroided spectrum.
    """
    spectrum = oms.MSSpectrum()
    spectrum.setMSLevel(1)
    spectrum.setRT(120.0)

    # Sparse peaks (centroid-like) with large gaps between them
    mz_values = [200.1, 300.2, 400.3, 450.4, 500.5, 550.6, 600.7, 700.8, 800.9, 900.0]
    intensity_values = [5000.0, 3000.0, 8000.0, 2000.0, 15000.0, 7000.0, 4000.0, 6000.0, 1000.0, 9000.0]

    spectrum.set_peaks([mz_values, intensity_values])
    return spectrum


@click.command(help="Estimate peak type (profile/centroid) for spectra in an mzML file.")
@click.option("--input", "input_path", required=True, help="Input mzML file")
@click.option("--output", "output_path", required=True, help="Output TSV report")
def main(input_path: str, output_path: str) -> None:
    result = estimate_peak_type(input_path, output_path)
    click.echo(
        f"Analyzed {result['total_spectra']} spectra: "
        + ", ".join(f"{k}={v}" for k, v in result['type_counts'].items())
    )
    click.echo(f"Report written to {output_path}")


if __name__ == "__main__":
    main()
