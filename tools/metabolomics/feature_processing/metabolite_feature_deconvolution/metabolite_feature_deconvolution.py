"""
Metabolite Feature Deconvolution
================================
Group metabolite features by adduct relationships (e.g. [M+H]+, [M+Na]+,
[M+K]+) into consensus feature groups. This helps identify features that
originate from the same underlying metabolite.

Wraps pyopenms.MetaboliteFeatureDeconvolution for adduct grouping.

Usage
-----
    python metabolite_feature_deconvolution.py --input features.featureXML --output grouped.consensusXML
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def deconvolve_adducts(
    input_path: str,
    output_path: str,
) -> int:
    """Group features by adduct relationships using MetaboliteFeatureDeconvolution.

    Parameters
    ----------
    input_path : str
        Path to input featureXML file.
    output_path : str
        Path to output consensusXML file with grouped features.

    Returns
    -------
    int
        Number of adduct groups (consensus features) found.
    """
    # Load input features
    in_fm = oms.FeatureMap()
    print(f"Loading {input_path} ...")
    oms.FeatureXMLFile().load(input_path, in_fm)

    # Run MetaboliteFeatureDeconvolution
    mfd = oms.MetaboliteFeatureDeconvolution()
    out_fm = oms.FeatureMap()
    out_cm = oms.ConsensusMap()
    out_cm_p = oms.ConsensusMap()

    mfd.compute(in_fm, out_fm, out_cm, out_cm_p)

    # Store the grouped consensus map
    oms.ConsensusXMLFile().store(output_path, out_cm)

    n_groups = out_cm.size()
    print(f"Found {n_groups} adduct groups -> {output_path}")
    return n_groups


def print_group_summary(output_path: str, top_n: int = 20) -> None:
    """Print a summary of the consensus groups.

    Parameters
    ----------
    output_path : str
        Path to the consensusXML file.
    top_n : int
        Number of groups to display.
    """
    cons_map = oms.ConsensusMap()
    oms.ConsensusXMLFile().load(output_path, cons_map)

    if cons_map.size() == 0:
        print("No adduct groups found.")
        return

    features = list(cons_map)
    features.sort(key=lambda f: f.getIntensity(), reverse=True)

    display = features[:top_n]
    print(
        f"\nTop {len(display)} adduct groups (by intensity):\n"
        f"{'#':>5}  {'RT (s)':>10}  {'m/z':>12}  {'Members':>8}  {'Intensity':>14}"
    )
    print("-" * 58)
    for i, cf in enumerate(display, 1):
        n_members = cf.size()
        print(
            f"{i:>5}  {cf.getRT():>10.2f}  {cf.getMZ():>12.4f}  "
            f"{n_members:>8}  {cf.getIntensity():>14.3e}"
        )


@click.command(help="Group metabolite features by adduct relationships.")
@click.option(
    "--input", "input_file", required=True, help="Input featureXML file"
)
@click.option("--output", required=True, help="Output consensusXML file")
@click.option(
    "--top",
    type=int,
    default=20,
    help="Number of top groups to print (default: 20)",
)
def main(input_file, output, top):
    n_groups = deconvolve_adducts(input_file, output)
    if n_groups > 0:
        print_group_summary(output, top)


if __name__ == "__main__":
    main()
