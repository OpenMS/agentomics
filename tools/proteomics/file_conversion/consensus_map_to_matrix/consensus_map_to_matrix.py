"""
Consensus Map to Matrix
=======================
Convert a consensusXML file to a quantification matrix (TSV).

Usage
-----
    python consensus_map_to_matrix.py --input consensus.consensusXML --output matrix.tsv
"""

import csv

import click
import pyopenms as oms


def load_consensus_map(input_path: str) -> oms.ConsensusMap:
    """Load a consensusXML file."""
    cmap = oms.ConsensusMap()
    oms.ConsensusXMLFile().load(input_path, cmap)
    return cmap


def consensus_to_matrix(input_path: str, output_path: str) -> dict:
    """Convert consensusXML to a quantification matrix TSV.

    Each row is a consensus feature. Columns include RT, MZ, charge,
    and intensity for each input map.

    Returns statistics about the conversion.
    """
    cmap = load_consensus_map(input_path)

    # Get column descriptions (file mappings)
    column_headers = cmap.getColumnHeaders()
    map_indices = sorted(column_headers.keys())
    map_labels = []
    for idx in map_indices:
        desc = column_headers[idx]
        label = desc.filename if desc.filename else f"map_{idx}"
        map_labels.append(label)

    n_maps = len(map_indices)
    rows = []

    for cf in cmap:
        row = {
            "rt": round(cf.getRT(), 4),
            "mz": round(cf.getMZ(), 6),
            "charge": cf.getCharge(),
            "quality": round(cf.getQuality(), 4),
        }

        # Initialize intensities for all maps
        intensities = {idx: 0.0 for idx in map_indices}
        for handle in cf.getFeatureList():
            map_idx = handle.getMapIndex()
            if map_idx in intensities:
                intensities[map_idx] = handle.getIntensity()

        for i, idx in enumerate(map_indices):
            row[f"intensity_{i}"] = round(intensities[idx], 4)

        rows.append(row)

    # Write TSV
    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        header = ["rt", "mz", "charge", "quality"] + [f"intensity_{i}" for i in range(n_maps)]
        writer.writerow(header)
        for row in rows:
            writer.writerow([row.get(col, "") for col in header])

    return {
        "consensus_features": len(rows),
        "n_maps": n_maps,
        "map_labels": map_labels,
    }


def create_synthetic_consensus(output_path: str, n_features: int = 5, n_maps: int = 3) -> None:
    """Create a synthetic consensusXML file for testing."""
    cmap = oms.ConsensusMap()

    # Set up column headers
    for i in range(n_maps):
        desc = oms.ColumnHeader()
        desc.filename = f"sample_{i}.mzML"
        desc.label = f"sample_{i}"
        desc.size = n_features
        cmap.setColumnHeaders({i: desc for i in range(n_maps)})

    # Re-set with proper dict
    headers = {}
    for i in range(n_maps):
        desc = oms.ColumnHeader()
        desc.filename = f"sample_{i}.mzML"
        desc.label = f"sample_{i}"
        desc.size = n_features
        headers[i] = desc
    cmap.setColumnHeaders(headers)

    for j in range(n_features):
        cf = oms.ConsensusFeature()
        cf.setRT(100.0 + j * 10)
        cf.setMZ(500.0 + j * 50)
        cf.setCharge(2)
        cf.setQuality(0.9)

        for i in range(n_maps):
            peak = oms.Peak2D()
            peak.setRT(100.0 + j * 10 + i * 0.1)
            peak.setMZ(500.0 + j * 50)
            peak.setIntensity(10000.0 + i * 1000 + j * 500)
            cf.insert(i, peak, j)
        cmap.push_back(cf)

    oms.ConsensusXMLFile().store(output_path, cmap)


@click.command(help="Convert consensusXML to quantification matrix.")
@click.option("--input", "input", required=True, help="Input consensusXML file")
@click.option("--output", required=True, help="Output TSV matrix file")
def main(input, output) -> None:
    stats = consensus_to_matrix(input, output)
    print(f"Exported {stats['consensus_features']} features across {stats['n_maps']} maps to {output}")


if __name__ == "__main__":
    main()
