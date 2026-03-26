"""
Multiplex Feature Finder
========================
Detect multiplex-labeled peptide features (e.g. SILAC light/heavy pairs)
from an mzML file using the FeatureFinderMultiplexAlgorithm.

Usage
-----
    python multiplex_feature_finder.py --input run.mzML --output features.featureXML \
        --labels "[][Lys8]"

Supported label names: Arg6, Arg10, Lys4, Lys6, Lys8, Leu3,
    Dimethyl0, Dimethyl4, Dimethyl6, Dimethyl8,
    ICPL0, ICPL4, ICPL6, ICPL10.

The labels string uses the OpenMS format: each channel is enclosed in brackets,
with comma-separated label names. The light (unlabeled) channel is [].
Examples:
    "[][Lys8]"              - SILAC light/heavy (Lys8)
    "[][Lys4][Lys8]"        - SILAC triple
    "[][Lys8,Arg10]"        - SILAC light/heavy (Lys8 + Arg10)
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

# Known labels and their mass shifts (Da)
KNOWN_LABELS = {
    "Arg6": 6.0201290268,
    "Arg10": 10.0082686,
    "Lys4": 4.0251069836,
    "Lys6": 6.0201290268,
    "Lys8": 8.0141988132,
    "Leu3": 3.01883,
    "Dimethyl0": 28.0313,
    "Dimethyl4": 32.056407,
    "Dimethyl6": 34.063117,
    "Dimethyl8": 36.07567,
    "ICPL0": 105.021464,
    "ICPL4": 109.046571,
    "ICPL6": 111.041593,
    "ICPL10": 115.0667,
}


def find_multiplex_features(
    input_path: str,
    output_path: str,
    labels: str = "[][Lys8]",
    charge_low: int = 1,
    charge_high: int = 4,
    mz_tolerance: float = 10.0,
    rt_typical: float = 40.0,
    intensity_cutoff: float = 1000.0,
) -> int:
    """Detect multiplex-labeled peptide features from an mzML file.

    Parameters
    ----------
    input_path : str
        Path to the input mzML file.
    output_path : str
        Path for the output featureXML file.
    labels : str
        Labels parameter in OpenMS format, e.g. "[][Lys8]" or "[][Lys8,Arg10]".
        Each channel is enclosed in brackets. The light channel is [].
    charge_low : int
        Minimum charge state to consider.
    charge_high : int
        Maximum charge state to consider.
    mz_tolerance : float
        m/z tolerance in ppm.
    rt_typical : float
        Typical elution profile width in seconds.
    intensity_cutoff : float
        Minimum intensity cutoff.

    Returns
    -------
    int
        Number of features detected.
    """
    # Load experiment
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    # Set up algorithm
    algo = oms.FeatureFinderMultiplexAlgorithm()
    params = algo.getDefaults()
    params.setValue("algorithm:labels", labels)
    params.setValue("algorithm:charge", f"{charge_low}:{charge_high}")
    params.setValue("algorithm:mz_tolerance", mz_tolerance)
    params.setValue("algorithm:mz_unit", "ppm")
    params.setValue("algorithm:rt_typical", rt_typical)
    params.setValue("algorithm:intensity_cutoff", intensity_cutoff)
    algo.setParameters(params)

    # Run feature detection
    algo.run(exp, False)

    # Get results
    feature_map = algo.getFeatureMap()

    # Store feature map
    oms.FeatureXMLFile().store(output_path, feature_map)

    return feature_map.size()


def create_synthetic_multiplex_mzml(
    output_path: str,
    light_mz: float = 500.0,
    mass_shift: float = 8.0142,
    charge: int = 2,
    rt_center: float = 100.0,
    intensity: float = 1e5,
    n_scans: int = 20,
) -> None:
    """Create a synthetic mzML with light/heavy peptide pairs.

    Parameters
    ----------
    output_path : str
        Path for the output mzML file.
    light_mz : float
        m/z of the light peptide.
    mass_shift : float
        Mass difference between light and heavy labels in Da.
    charge : int
        Charge state of the peptides.
    rt_center : float
        Center of the elution profile in seconds.
    intensity : float
        Peak intensity.
    n_scans : int
        Number of MS1 scans to generate.
    """
    exp = oms.MSExperiment()
    heavy_mz = light_mz + mass_shift / charge

    for i in range(n_scans):
        rt = rt_center - n_scans / 2 + i
        spec = oms.MSSpectrum()
        spec.setRT(rt)
        spec.setMSLevel(1)

        # Gaussian-like intensity profile
        dist = abs(rt - rt_center)
        scale = max(0.01, 1.0 - (dist / (n_scans / 2)) ** 2)
        inten = intensity * scale

        # Add isotope envelope for light peptide
        for iso in range(4):
            mz_iso = light_mz + iso * 1.003355 / charge
            peak = oms.Peak1D()
            peak.setMZ(mz_iso)
            peak.setIntensity(inten * (0.6 ** iso))
            spec.push_back(peak)

        # Add isotope envelope for heavy peptide
        for iso in range(4):
            mz_iso = heavy_mz + iso * 1.003355 / charge
            peak = oms.Peak1D()
            peak.setMZ(mz_iso)
            peak.setIntensity(inten * (0.6 ** iso))
            spec.push_back(peak)

        spec.sortByPosition()
        exp.addSpectrum(spec)

    exp.sortSpectra()
    oms.MzMLFile().store(output_path, exp)


@click.command(help="Detect multiplex-labeled peptide features from an mzML file.")
@click.option("--input", "input_path", required=True, help="Input mzML file")
@click.option("--output", "output_path", required=True, help="Output featureXML file")
@click.option(
    "--labels",
    default="[][Lys8]",
    help="Labels in OpenMS format, e.g. '[][Lys8]' or '[][Lys8,Arg10]' (default: [][Lys8])",
)
@click.option("--charge-low", default=1, help="Minimum charge state (default: 1)")
@click.option("--charge-high", default=4, help="Maximum charge state (default: 4)")
@click.option("--mz-tolerance", default=10.0, help="m/z tolerance in ppm (default: 10.0)")
@click.option("--rt-typical", default=40.0, help="Typical RT width in seconds (default: 40.0)")
@click.option(
    "--intensity-cutoff", default=1000.0, help="Intensity cutoff (default: 1000.0)"
)
def main(
    input_path,
    output_path,
    labels,
    charge_low,
    charge_high,
    mz_tolerance,
    rt_typical,
    intensity_cutoff,
) -> None:
    n_features = find_multiplex_features(
        input_path,
        output_path,
        labels=labels,
        charge_low=charge_low,
        charge_high=charge_high,
        mz_tolerance=mz_tolerance,
        rt_typical=rt_typical,
        intensity_cutoff=intensity_cutoff,
    )
    click.echo(f"Detected {n_features} multiplex features, saved to {output_path}")


if __name__ == "__main__":
    main()
