"""
Feature Linker Labeled
======================
Link labeled feature pairs (e.g., SILAC light/heavy) within a single
LC-MS/MS run using the labeled grouping algorithm.

Reads a single featureXML file and produces a consensusXML file where
paired features (e.g., light and heavy isotope labels) are grouped.

Usage
-----
    python feature_linker_labeled.py --input run.featureXML \
        --output consensus.consensusXML --rt-estimate 1000.0 --mz-estimate 5.0
"""

import sys
from typing import List, Optional

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def link_labeled_features(
    input_path: str,
    output_path: str,
    rt_estimate: float = 1000.0,
    mz_estimate: float = 5.0,
    mz_pair_dists: Optional[List[float]] = None,
    rt_pair_dist: float = -20.0,
) -> int:
    """Link labeled feature pairs within a single run.

    Parameters
    ----------
    input_path : str
        Path to the input featureXML file.
    output_path : str
        Path for the output consensusXML file.
    rt_estimate : float
        Maximum allowed RT deviation for pairing (default 1000.0).
        Sets both rt_dev_low and rt_dev_high internally.
    mz_estimate : float
        Maximum allowed m/z deviation around the expected pair distance
        (default 5.0).
    mz_pair_dists : list of float, optional
        Expected m/z distances between labeled pairs.  When *None* the
        pyopenms default ([4.0]) is used.  For SILAC K+8 pass [8.0142].
    rt_pair_dist : float
        Expected RT distance between paired features (default -20.0).
        Negative values indicate the heavier partner elutes earlier.

    Returns
    -------
    int
        Number of consensus features (linked pairs) in the result.
    """
    # Load the feature map
    fm = oms.FeatureMap()
    oms.FeatureXMLFile().load(input_path, fm)

    # Ensure all features have unique IDs (required by the algorithm)
    fm.ensureUniqueId()
    for i in range(fm.size()):
        fm[i].ensureUniqueId()

    # Configure the linker
    linker = oms.FeatureGroupingAlgorithmLabeled()
    params = linker.getDefaults()
    params.setValue("rt_estimate", "false")
    params.setValue("rt_pair_dist", rt_pair_dist)
    params.setValue("rt_dev_low", rt_estimate)
    params.setValue("rt_dev_high", rt_estimate)
    params.setValue("mz_dev", mz_estimate)
    if mz_pair_dists is not None:
        params.setValue("mz_pair_dists", mz_pair_dists)
    linker.setParameters(params)

    # Set up consensus map with two column headers (light and heavy channels)
    consensus = oms.ConsensusMap()
    file_descriptions = consensus.getColumnHeaders()
    desc_light = oms.ColumnHeader()
    desc_light.filename = input_path
    desc_light.size = fm.size()
    desc_light.label = "light"
    file_descriptions[0] = desc_light
    desc_heavy = oms.ColumnHeader()
    desc_heavy.filename = input_path
    desc_heavy.size = fm.size()
    desc_heavy.label = "heavy"
    file_descriptions[1] = desc_heavy
    consensus.setColumnHeaders(file_descriptions)

    # Run the grouping (labeled takes a single map in a list)
    linker.group([fm], consensus)

    # Store result
    oms.ConsensusXMLFile().store(output_path, consensus)

    return consensus.size()


@click.command()
@click.option(
    "--input",
    "input_path",
    required=True,
    help="Input featureXML file.",
)
@click.option(
    "--output",
    required=True,
    help="Output consensusXML file path.",
)
@click.option(
    "--rt-estimate",
    default=1000.0,
    show_default=True,
    help="Maximum allowed RT deviation for pairing.",
)
@click.option(
    "--mz-estimate",
    default=5.0,
    show_default=True,
    help="Maximum allowed m/z deviation around expected pair distance.",
)
def main(input_path, output, rt_estimate, mz_estimate):
    """Link labeled feature pairs within a single LC-MS/MS run."""
    count = link_labeled_features(
        input_path, output, rt_estimate=rt_estimate, mz_estimate=mz_estimate
    )
    click.echo(f"Linked features into {count} consensus features.")
    click.echo(f"Output written to {output}")


if __name__ == "__main__":
    main()
