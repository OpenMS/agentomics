"""
Elution Peak Detector
=====================
Detect elution peaks from LC-MS data by first finding mass traces and then
splitting them into individual chromatographic peaks.

Wraps pyopenms.MassTraceDetection + pyopenms.ElutionPeakDetection to perform
a two-step workflow: mass trace detection followed by elution peak detection.

Usage
-----
    python elution_peak_detector.py --input run.mzML --output peaks.mzML \
        --width-filtering auto
"""

import math
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def detect_elution_peaks(
    input_path: str,
    output_path: str,
    width_filtering: str = "auto",
    mass_error_ppm: float = 10.0,
    noise_threshold: float = 1000.0,
) -> int:
    """Detect elution peaks in an LC-MS mzML file.

    Runs MassTraceDetection internally to find mass traces, then
    ElutionPeakDetection to split them into individual peaks.

    Parameters
    ----------
    input_path : str
        Path to input mzML file with MS1 spectra.
    output_path : str
        Path to output mzML file (peaks stored as chromatograms).
    width_filtering : str
        Width filtering mode: 'auto', 'fixed', or 'off' (default 'auto').
    mass_error_ppm : float
        Mass error tolerance in ppm for mass trace detection (default 10.0).
    noise_threshold : float
        Noise threshold for mass trace detection (default 1000.0).

    Returns
    -------
    int
        Number of elution peaks detected.
    """
    # Load experiment
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    # Step 1: detect mass traces
    mtd = oms.MassTraceDetection()
    mtd_params = mtd.getDefaults()
    mtd_params.setValue("mass_error_ppm", mass_error_ppm)
    mtd_params.setValue("noise_threshold_int", noise_threshold)
    mtd.setParameters(mtd_params)

    mass_traces = []
    mtd.run(exp, mass_traces, 0)

    if not mass_traces:
        # No mass traces found, write empty output
        out_exp = oms.MSExperiment()
        oms.MzMLFile().store(output_path, out_exp)
        return 0

    # Step 2: detect elution peaks
    epd = oms.ElutionPeakDetection()
    epd_params = epd.getDefaults()
    epd_params.setValue("width_filtering", width_filtering)
    epd.setParameters(epd_params)

    split_traces = []
    epd.detectPeaks(mass_traces, split_traces)

    # Store elution peaks as chromatograms in output mzML
    out_exp = oms.MSExperiment()
    for i, mt in enumerate(split_traces):
        chrom = oms.MSChromatogram()
        chrom.setNativeID(f"elution_peak_{i}")

        # Set product m/z from the trace centroid
        product = oms.Product()
        product.setMZ(mt.getCentroidMZ())
        chrom.setProduct(product)

        # Extract RT/intensity from convex hull
        hull = mt.getConvexhull()
        points = hull.getHullPoints()

        rt_set = set()
        for pt in points:
            rt_set.add(pt[0])

        sorted_rts = sorted(rt_set)
        if sorted_rts:
            rts = []
            intensities = []
            for rt in sorted_rts:
                rts.append(rt)
                intensities.append(float(mt.getIntensity(False)))
            chrom.set_peaks((rts, intensities))

        out_exp.addChromatogram(chrom)

    oms.MzMLFile().store(output_path, out_exp)

    return len(split_traces)


def create_synthetic_lcms_mzml(output_path: str, n_scans: int = 20) -> None:
    """Create a synthetic LC-MS mzML with a Gaussian elution profile at m/z=500.

    The peak has a bell-curve intensity centered at scan n_scans//2.
    """
    exp = oms.MSExperiment()
    for i in range(n_scans):
        spec = oms.MSSpectrum()
        spec.setMSLevel(1)
        spec.setRT(float(i) * 5.0)
        # Gaussian intensity profile
        intensity = 10000.0 * math.exp(-0.5 * ((i - n_scans // 2) / 3.0) ** 2)
        spec.set_peaks(([500.0], [intensity]))
        exp.addSpectrum(spec)
    oms.MzMLFile().store(output_path, exp)


@click.command(help="Detect elution peaks in LC-MS data.")
@click.option("--input", "input_path", required=True, help="Input mzML file")
@click.option("--output", "output_path", required=True, help="Output mzML file with peaks as chromatograms")
@click.option(
    "--width-filtering",
    default="auto",
    type=click.Choice(["auto", "fixed", "off"]),
    help="Peak width filtering mode",
)
@click.option("--mass-error-ppm", default=10.0, type=float, help="Mass error in ppm")
@click.option("--noise-threshold", default=1000.0, type=float, help="Noise intensity threshold")
def main(input_path, output_path, width_filtering, mass_error_ppm, noise_threshold) -> None:
    peaks = detect_elution_peaks(
        input_path, output_path, width_filtering, mass_error_ppm, noise_threshold
    )
    click.echo(f"Detected {peaks} elution peaks, saved to {output_path}")


if __name__ == "__main__":
    main()
