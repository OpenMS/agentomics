"""
Mass Trace Detector
===================
Detect mass traces in LC-MS data using MassTraceDetection.

A mass trace is a series of connected peaks across consecutive spectra
sharing a similar m/z value. This tool wraps pyopenms.MassTraceDetection
and stores detected traces as chromatograms in an output mzML file.

Usage
-----
    python mass_trace_detector.py --input run.mzML --output traces.mzML \
        --mass-error-ppm 10.0 --noise-threshold 1000
"""

import math
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def detect_mass_traces(
    input_path: str,
    output_path: str,
    mass_error_ppm: float = 10.0,
    noise_threshold: float = 1000.0,
) -> int:
    """Detect mass traces in an LC-MS mzML file.

    Parameters
    ----------
    input_path : str
        Path to input mzML file with MS1 spectra.
    output_path : str
        Path to output mzML file (traces stored as chromatograms).
    mass_error_ppm : float
        Mass error tolerance in ppm (default 10.0).
    noise_threshold : float
        Noise threshold for intensity filtering (default 1000.0).

    Returns
    -------
    int
        Number of mass traces detected.
    """
    # Load experiment
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    # Configure mass trace detection
    mtd = oms.MassTraceDetection()
    params = mtd.getDefaults()
    params.setValue("mass_error_ppm", mass_error_ppm)
    params.setValue("noise_threshold_int", noise_threshold)
    mtd.setParameters(params)

    # Run detection
    mass_traces = []
    mtd.run(exp, mass_traces, 0)

    # Store mass trace data as chromatograms in output mzML
    out_exp = oms.MSExperiment()
    for i, mt in enumerate(mass_traces):
        chrom = oms.MSChromatogram()
        chrom.setNativeID(f"mass_trace_{i}")

        # Extract RT/intensity data from the convex hull
        hull = mt.getConvexhull()
        points = hull.getHullPoints()

        # Group points by RT to extract the trace profile
        rt_intensity = {}
        for pt in points:
            rt = pt[0]
            mz = pt[1]
            # Use the RT as key; convex hull points form a boundary
            if rt not in rt_intensity:
                rt_intensity[rt] = mz
            else:
                # Keep the higher value as intensity proxy
                rt_intensity[rt] = max(rt_intensity[rt], mz)

        # Set chromatogram peaks from sorted RT values
        sorted_rts = sorted(rt_intensity.keys())
        if sorted_rts:
            # Store the centroid mz as product m/z
            product = oms.Product()
            product.setMZ(mt.getCentroidMZ())
            chrom.setProduct(product)

            rts = []
            intensities = []
            for rt in sorted_rts:
                rts.append(rt)
                intensities.append(float(mt.getIntensity(False)))
            chrom.set_peaks((rts, intensities))

        out_exp.addChromatogram(chrom)

    oms.MzMLFile().store(output_path, out_exp)

    return len(mass_traces)


def create_synthetic_lcms_mzml(output_path: str, n_scans: int = 20) -> None:
    """Create a synthetic LC-MS mzML with a persistent peak at m/z=500.

    The peak has a Gaussian intensity profile centered at scan 10.
    """
    exp = oms.MSExperiment()
    for i in range(n_scans):
        spec = oms.MSSpectrum()
        spec.setMSLevel(1)
        spec.setRT(float(i) * 5.0)
        intensity = 10000.0 * math.exp(-0.5 * ((i - n_scans // 2) / 3.0) ** 2)
        spec.set_peaks(([500.0], [intensity]))
        exp.addSpectrum(spec)
    oms.MzMLFile().store(output_path, exp)


@click.command(help="Detect mass traces in LC-MS data.")
@click.option("--input", "input_path", required=True, help="Input mzML file")
@click.option("--output", "output_path", required=True, help="Output mzML file with traces as chromatograms")
@click.option("--mass-error-ppm", default=10.0, type=float, help="Mass error in ppm")
@click.option("--noise-threshold", default=1000.0, type=float, help="Noise intensity threshold")
def main(input_path, output_path, mass_error_ppm, noise_threshold) -> None:
    traces = detect_mass_traces(input_path, output_path, mass_error_ppm, noise_threshold)
    click.echo(f"Detected {traces} mass traces, saved to {output_path}")


if __name__ == "__main__":
    main()
