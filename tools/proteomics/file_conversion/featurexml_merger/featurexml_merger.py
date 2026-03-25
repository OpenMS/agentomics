"""
featureXML Merger
=================
Merge multiple featureXML files into a single featureXML file.

Usage
-----
    python featurexml_merger.py --inputs f1.featureXML f2.featureXML --output merged.featureXML
"""

from typing import List

import click
import pyopenms as oms


def load_featurexml(input_path: str) -> oms.FeatureMap:
    """Load a featureXML file."""
    fm = oms.FeatureMap()
    oms.FeatureXMLFile().load(input_path, fm)
    return fm


def save_featurexml(feature_map: oms.FeatureMap, output_path: str) -> None:
    """Save a FeatureMap to featureXML."""
    oms.FeatureXMLFile().store(output_path, feature_map)


def merge_feature_maps(input_paths: List[str], output_path: str) -> dict:
    """Merge multiple featureXML files into one.

    Returns statistics about the merge.
    """
    merged = oms.FeatureMap()
    file_counts = {}

    for path in input_paths:
        fm = load_featurexml(path)
        count = fm.size()
        file_counts[path] = count
        for feature in fm:
            merged.push_back(feature)

    # Sort by RT
    merged.sortByRT()

    save_featurexml(merged, output_path)

    return {
        "file_counts": file_counts,
        "total_features": merged.size(),
    }


def create_synthetic_featurexml(output_path: str, n_features: int = 5, rt_offset: float = 0.0) -> None:
    """Create a synthetic featureXML file for testing."""
    fm = oms.FeatureMap()
    for i in range(n_features):
        f = oms.Feature()
        f.setRT(100.0 + rt_offset + i * 10)
        f.setMZ(500.0 + i * 50)
        f.setIntensity(10000.0 + i * 1000)
        f.setCharge(2)
        f.setOverallQuality(0.9)
        fm.push_back(f)
    save_featurexml(fm, output_path)


@click.command(help="Merge multiple featureXML files.")
@click.option("--inputs", multiple=True, required=True, help="Input featureXML files")
@click.option("--output", required=True, help="Output merged featureXML file")
def main(inputs, output) -> None:
    stats = merge_feature_maps(list(inputs), output)
    print(f"Merged {stats['total_features']} features from {len(stats['file_counts'])} files to {output}")


if __name__ == "__main__":
    main()
