"""
RT Alignment - Pose Clustering
===============================
Align retention times between feature maps using the pose clustering algorithm.

This tool aligns an input featureXML file to a reference featureXML file by
computing a retention time transformation based on pose clustering of matching
features.

Usage
-----
    python rt_alignment_pose_clustering.py \\
        --reference ref.featureXML --input sample.featureXML --output aligned.featureXML
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def create_synthetic_featurexml(
    output_path: str, n_features: int = 10, rt_offset: float = 0.0
) -> None:
    """Create a synthetic featureXML file for testing.

    Parameters
    ----------
    output_path : str
        Path to write the featureXML file.
    n_features : int
        Number of features to generate.
    rt_offset : float
        Offset to add to retention times (seconds).
    """
    fm = oms.FeatureMap()
    fm.setUniqueId(1)
    for i in range(n_features):
        f = oms.Feature()
        f.setMZ(500.0 + i * 100.0)
        f.setRT(100.0 + i * 100.0 + rt_offset)
        f.setIntensity(1e5)
        f.setCharge(2)
        f.setOverallQuality(1.0)
        f.setUniqueId(i + 2)
        fm.push_back(f)
    oms.FeatureXMLFile().store(output_path, fm)


def align_pose_clustering(
    reference_path: str, input_path: str, output_path: str
) -> dict:
    """Align an input feature map to a reference using pose clustering.

    Parameters
    ----------
    reference_path : str
        Path to the reference featureXML file.
    input_path : str
        Path to the input featureXML file to be aligned.
    output_path : str
        Path to write the aligned featureXML file.

    Returns
    -------
    dict
        Statistics about the alignment:
        - n_features: number of features in the input map
        - max_rt_shift: maximum RT shift applied (seconds)
        - n_trafo_points: number of transformation data points
    """
    # Load reference and input feature maps
    reference_map = oms.FeatureMap()
    oms.FeatureXMLFile().load(reference_path, reference_map)

    input_map = oms.FeatureMap()
    oms.FeatureXMLFile().load(input_path, input_map)

    # Record original RTs for statistics
    original_rts = [input_map[i].getRT() for i in range(input_map.size())]

    # Set up and run the pose clustering alignment
    aligner = oms.MapAlignmentAlgorithmPoseClustering()
    aligner.setReference(reference_map)

    trafo = oms.TransformationDescription()
    aligner.align(input_map, trafo)

    # Apply the transformation to the input map
    transformer = oms.MapAlignmentTransformer()
    transformer.transformRetentionTimes(input_map, trafo, False)

    # Compute max RT shift
    aligned_rts = [input_map[i].getRT() for i in range(input_map.size())]
    max_rt_shift = 0.0
    for orig, aligned in zip(original_rts, aligned_rts):
        shift = abs(orig - aligned)
        if shift > max_rt_shift:
            max_rt_shift = shift

    # Save the aligned map
    oms.FeatureXMLFile().store(output_path, input_map)

    return {
        "n_features": input_map.size(),
        "max_rt_shift": max_rt_shift,
        "n_trafo_points": len(trafo.getDataPoints()),
    }


@click.command(help="Align retention times using pose clustering.")
@click.option(
    "--reference", required=True, help="Reference featureXML file."
)
@click.option(
    "--input", "input_path", required=True, help="Input featureXML file to align."
)
@click.option(
    "--output", required=True, help="Output aligned featureXML file."
)
def main(reference: str, input_path: str, output: str) -> None:
    """CLI entry point for pose clustering alignment."""
    stats = align_pose_clustering(reference, input_path, output)
    click.echo(
        f"Aligned {stats['n_features']} features "
        f"(max RT shift: {stats['max_rt_shift']:.2f}s, "
        f"{stats['n_trafo_points']} transformation points) -> {output}"
    )


if __name__ == "__main__":
    main()
