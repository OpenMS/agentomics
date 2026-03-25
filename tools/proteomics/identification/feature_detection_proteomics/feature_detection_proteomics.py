"""
Feature Detection for Proteomics LC-MS Data
============================================
Detect peptide features (isotope envelopes) in an mzML file using the
pyopenms FeatureFinderCentroided algorithm.  Results are written to a
featureXML file which can be opened in TOPPView.

Usage
-----
    python feature_detection_proteomics.py --input sample.mzML
    python feature_detection_proteomics.py --input sample.mzML --output features.featureXML
"""


import click
import pyopenms as oms


def detect_features(
    input_path: str,
    output_path: str,
) -> oms.FeatureMap:
    """Run FeatureFinderCentroided on an mzML file.

    Parameters
    ----------
    input_path:
        Path to the centroided mzML file.
    output_path:
        Path for the output featureXML file.

    Returns
    -------
    pyopenms.FeatureMap
        Map of detected features.
    """
    exp = oms.MSExperiment()
    print(f"Loading {input_path} …")
    oms.MzMLFile().load(input_path, exp)
    exp.updateRanges()

    feature_map = oms.FeatureMap()
    seeds = oms.FeatureMap()

    ff = oms.FeatureFinderAlgorithmPicked()
    params = ff.getParameters()
    ff.run(exp, feature_map, params, seeds)

    feature_map.setUniqueIds()
    oms.FeatureXMLFile().store(output_path, feature_map)
    print(f"Detected {feature_map.size()} features → {output_path}")
    return feature_map


def print_feature_summary(feature_map: oms.FeatureMap) -> None:
    """Print a tabular summary of detected features."""
    if feature_map.size() == 0:
        print("No features detected.")
        return

    print(f"\n{'#':>5}  {'RT (s)':>10}  {'m/z':>12}  {'Charge':>6}  {'Intensity':>14}")
    print("-" * 56)
    for i, feature in enumerate(feature_map, 1):
        print(
            f"{i:>5}  {feature.getRT():>10.2f}  {feature.getMZ():>12.4f}  "
            f"{feature.getCharge():>6}  {feature.getIntensity():>14.3e}"
        )


@click.command(help="Detect peptide features in an mzML file using pyopenms.")
@click.option("--input", "input", required=True, help="Centroided mzML input file")
@click.option("--output", default=None, help="Output featureXML file (default: <input>.featureXML)")
def main(input, output):
    output_path = output or input.replace(".mzML", ".featureXML")
    feature_map = detect_features(input, output_path)
    print_feature_summary(feature_map)


if __name__ == "__main__":
    main()
