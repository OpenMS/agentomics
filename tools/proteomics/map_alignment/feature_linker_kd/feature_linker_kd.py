"""
Feature Linker KD
=================
Link features across multiple LC-MS/MS runs using the KD-tree-based
grouping algorithm.

Reads multiple featureXML files and produces a single consensusXML file
where matching features (similar m/z and RT) are grouped together.

Usage
-----
    python feature_linker_kd.py --inputs run1.featureXML run2.featureXML \
        --output consensus.consensusXML --rt-tol 30.0 --mz-tol 0.01
"""

import sys
from typing import List

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def link_features_kd(
    input_paths: List[str],
    output_path: str,
    rt_tol: float = 30.0,
    mz_tol: float = 0.01,
) -> int:
    """Link features across multiple runs using the KD-tree algorithm.

    Parameters
    ----------
    input_paths : List[str]
        Paths to featureXML files to link.
    output_path : str
        Path for the output consensusXML file.
    rt_tol : float
        Maximum allowed RT difference in seconds (default 30.0).
    mz_tol : float
        Maximum allowed m/z difference in Da (default 0.01).

    Returns
    -------
    int
        Number of consensus features in the result.
    """
    # Load feature maps
    feature_maps = []
    for path in input_paths:
        fm = oms.FeatureMap()
        oms.FeatureXMLFile().load(path, fm)
        feature_maps.append(fm)

    # Configure the linker
    linker = oms.FeatureGroupingAlgorithmKD()
    params = linker.getDefaults()
    params.setValue("warp:enabled", "false")
    params.setValue("link:rt_tol", rt_tol)
    params.setValue("link:mz_tol", mz_tol)
    linker.setParameters(params)

    # Set up consensus map with column headers
    consensus = oms.ConsensusMap()
    file_descriptions = consensus.getColumnHeaders()
    for i, fm in enumerate(feature_maps):
        desc = oms.ColumnHeader()
        desc.filename = input_paths[i]
        desc.size = fm.size()
        file_descriptions[i] = desc
    consensus.setColumnHeaders(file_descriptions)

    # Run the grouping
    linker.group(feature_maps, consensus)

    # Store result
    oms.ConsensusXMLFile().store(output_path, consensus)

    return consensus.size()


@click.command()
@click.option(
    "--inputs",
    required=True,
    multiple=True,
    help="Input featureXML files (specify multiple times or space-separated).",
)
@click.option(
    "--output",
    required=True,
    help="Output consensusXML file path.",
)
@click.option(
    "--rt-tol",
    default=30.0,
    show_default=True,
    help="Maximum RT difference in seconds.",
)
@click.option(
    "--mz-tol",
    default=0.01,
    show_default=True,
    help="Maximum m/z difference in Da.",
)
def main(inputs, output, rt_tol, mz_tol):
    """Link features across multiple LC-MS/MS runs using KD-tree grouping."""
    input_paths = list(inputs)
    if len(input_paths) < 2:
        click.echo("Error: At least two input files are required.", err=True)
        sys.exit(1)

    count = link_features_kd(input_paths, output, rt_tol=rt_tol, mz_tol=mz_tol)
    click.echo(f"Linked features into {count} consensus features.")
    click.echo(f"Output written to {output}")


if __name__ == "__main__":
    main()
