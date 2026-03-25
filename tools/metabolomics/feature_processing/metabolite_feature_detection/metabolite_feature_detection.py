"""
Metabolite Feature Detection
=============================
Detect small-molecule features (isotope envelopes) in an LC-MS mzML file
using the pyopenms FeatureFinderMetabo algorithm.  Results are written to a
featureXML file which can be opened in TOPPView.

Usage
-----
    python metabolite_feature_detection.py --input sample.mzML
    python metabolite_feature_detection.py --input sample.mzML --output features.featureXML --noise 1e5
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit(
        "pyopenms is required. Install it with:  pip install pyopenms"
    )


def detect_metabolite_features(
    input_path: str,
    output_path: str,
    noise_threshold: float = 1e4,
) -> oms.FeatureMap:
    """Run FeatureFinderMetabo on an mzML file.

    Parameters
    ----------
    input_path:
        Path to the centroided mzML file.
    output_path:
        Path for the output featureXML file.
    noise_threshold:
        Minimum peak intensity to consider during mass tracing (default 1e4).

    Returns
    -------
    pyopenms.FeatureMap
        Map of detected metabolite features.
    """
    exp = oms.MSExperiment()
    print(f"Loading {input_path} …")
    oms.MzMLFile().load(input_path, exp)
    exp.updateRanges()

    # --- Mass tracing ---
    mass_traces = []
    mt_params = oms.MassTraceDetection().getDefaults()
    mt_params.setValue("noise_threshold_int", noise_threshold)
    mt_det = oms.MassTraceDetection()
    mt_det.setParameters(mt_params)
    mt_det.run(exp, mass_traces, 0)
    print(f"Mass traces found: {len(mass_traces)}")

    # --- Elution peak detection ---
    mass_traces_split = []
    mass_traces_final = []
    epd_params = oms.ElutionPeakDetection().getDefaults()
    epd = oms.ElutionPeakDetection()
    epd.setParameters(epd_params)
    epd.detectPeaks(mass_traces, mass_traces_split)
    if epd.getParameters().getValue("width_filtering") == "auto":
        epd.filterByPeakWidth(mass_traces_split, mass_traces_final)
    else:
        mass_traces_final = mass_traces_split

    # --- Feature detection ---
    feature_map = oms.FeatureMap()
    output_chromatograms = []  # pyopenms >= 3.x requires a list here
    ffm_params = oms.FeatureFindingMetabo().getDefaults()
    ffm = oms.FeatureFindingMetabo()
    ffm.setParameters(ffm_params)
    ffm.run(mass_traces_final, feature_map, output_chromatograms)

    feature_map.setUniqueIds()
    oms.FeatureXMLFile().store(output_path, feature_map)
    print(f"Detected {feature_map.size()} metabolite features → {output_path}")
    return feature_map


def print_feature_summary(feature_map: oms.FeatureMap, top_n: int = 20) -> None:
    """Print a tabular summary of the top-N most intense features."""
    if feature_map.size() == 0:
        print("No features detected.")
        return

    features = list(feature_map)
    features.sort(key=lambda f: f.getIntensity(), reverse=True)

    display = features[:top_n]
    print(
        f"\nTop {len(display)} features (by intensity):\n"
        f"{'#':>5}  {'RT (s)':>10}  {'m/z':>12}  {'Charge':>6}  {'Intensity':>14}"
    )
    print("-" * 56)
    for i, feature in enumerate(display, 1):
        print(
            f"{i:>5}  {feature.getRT():>10.2f}  {feature.getMZ():>12.4f}  "
            f"{feature.getCharge():>6}  {feature.getIntensity():>14.3e}"
        )


@click.command()
@click.option("--input", "input_file", required=True, help="Centroided mzML input file")
@click.option("--output", default=None, help="Output featureXML file (default: <input>.featureXML)")
@click.option("--noise", type=float, default=1e4, help="Noise intensity threshold for mass tracing (default: 1e4)")
@click.option("--top", type=int, default=20, help="Number of top features to print (default: 20)")
def main(input_file, output, noise, top):
    output_path = output or input_file.replace(".mzML", "_metabolites.featureXML")
    feature_map = detect_metabolite_features(input_file, output_path, noise)
    print_feature_summary(feature_map, top)


if __name__ == "__main__":
    main()
