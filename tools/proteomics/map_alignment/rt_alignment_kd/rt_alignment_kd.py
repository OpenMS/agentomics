"""
RT Alignment - KD Tree
======================
Align retention times across multiple feature maps using a kd-tree based
reference-free algorithm.

This tool uses the KD-tree alignment algorithm to efficiently find
corresponding features across multiple LC-MS runs and compute LOWESS-based
RT transformations that minimize systematic RT shifts between runs.

Usage
-----
    python rt_alignment_kd.py --inputs run1.featureXML run2.featureXML --output-dir aligned/
    python rt_alignment_kd.py --inputs run1.featureXML run2.featureXML --output-dir aligned/ --model lowess
"""

import os
import sys
from collections import defaultdict
from typing import List

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def create_synthetic_featurexml(
    output_path: str, n_features: int = 30, rt_offset: float = 0.0
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
        f.setMZ(500.0 + i * 10.0)
        f.setRT(100.0 + i * 30.0 + rt_offset)
        f.setIntensity(1e5)
        f.setCharge(2)
        f.setOverallQuality(1.0)
        f.setUniqueId(i + 2)
        fm.push_back(f)
    oms.FeatureXMLFile().store(output_path, fm)


def align_kd(
    input_paths: List[str],
    output_dir: str,
    model: str = "lowess",
) -> dict:
    """Align multiple feature maps using the KD-tree algorithm.

    Uses the KD-tree based alignment from FeatureGroupingAlgorithmKD to find
    corresponding features across runs without requiring a designated reference.
    Computes RT transformations from the correspondences and applies them to
    each input map.

    Parameters
    ----------
    input_paths : list of str
        Paths to input featureXML files (at least 2).
    output_dir : str
        Directory to write the aligned featureXML files.
    model : str
        Transformation model: "lowess" or "linear".

    Returns
    -------
    dict
        Statistics about the alignment:
        - n_maps: number of input maps
        - n_features_per_map: list of feature counts per map
        - output_dir: path to the output directory
    """
    os.makedirs(output_dir, exist_ok=True)

    # Load all feature maps
    feature_maps = []
    n_features_per_map = []
    for path in input_paths:
        fm = oms.FeatureMap()
        oms.FeatureXMLFile().load(path, fm)
        feature_maps.append(fm)
        n_features_per_map.append(fm.size())

    # Run KD-tree feature grouping with warping enabled
    fga = oms.FeatureGroupingAlgorithmKD()
    params = fga.getDefaults()
    params.setValue(b"warp:enabled", "true")
    params.setValue(b"warp:rt_tol", 100.0)
    params.setValue(b"warp:mz_tol", 10.0)
    params.setValue(b"warp:min_rel_cc_size", 0.1)
    fga.setParameters(params)

    consensus = oms.ConsensusMap()
    fga.group(feature_maps, consensus)

    # Derive per-map transformations from the consensus map
    # Each consensus feature contains handles with original RTs from each map;
    # the consensus RT is the aligned (average) RT.
    trafo_data = defaultdict(list)
    for cf in consensus:
        consensus_rt = cf.getRT()
        for handle in cf.getFeatureList():
            map_idx = handle.getMapIndex()
            original_rt = handle.getRT()
            trafo_data[map_idx].append((original_rt, consensus_rt))

    # Build and apply transformations
    transformer = oms.MapAlignmentTransformer()
    for map_idx in range(len(feature_maps)):
        if map_idx in trafo_data and trafo_data[map_idx]:
            data_points = sorted(trafo_data[map_idx])
            trafo = oms.TransformationDescription()
            dp_list = [
                oms.TM_DataPoint(orig, target, "")
                for orig, target in data_points
            ]
            trafo.setDataPoints(dp_list)

            fit_model = b"lowess" if model == "lowess" else model.encode()
            trafo.fitModel(fit_model)
            transformer.transformRetentionTimes(
                feature_maps[map_idx], trafo, False
            )

    # Save aligned maps
    for path, fm in zip(input_paths, feature_maps):
        basename = os.path.basename(path)
        out_path = os.path.join(output_dir, basename)
        oms.FeatureXMLFile().store(out_path, fm)

    return {
        "n_maps": len(feature_maps),
        "n_features_per_map": n_features_per_map,
        "output_dir": output_dir,
    }


@click.command(help="Align retention times across feature maps using KD-tree algorithm.")
@click.option(
    "--inputs",
    multiple=True,
    required=True,
    help="Input featureXML files (at least 2).",
)
@click.option(
    "--output-dir",
    required=True,
    help="Directory for aligned featureXML output files.",
)
@click.option(
    "--model",
    default="lowess",
    type=click.Choice(["lowess", "linear"]),
    help="Transformation model (default: lowess).",
)
def main(inputs: tuple, output_dir: str, model: str) -> None:
    """CLI entry point for KD-tree based RT alignment."""
    stats = align_kd(list(inputs), output_dir, model=model)
    click.echo(
        f"Aligned {stats['n_maps']} maps "
        f"(features per map: {stats['n_features_per_map']}) -> {stats['output_dir']}"
    )


if __name__ == "__main__":
    main()
