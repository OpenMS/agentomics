"""
Consensus Map Normalizer
========================
Normalize intensity values across channels in a consensusXML file using
median or quantile normalization.

Usage
-----
    python consensus_map_normalizer.py --input consensus.consensusXML \
        --output normalized.consensusXML --method median
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

VALID_METHODS = ("median", "quantile")


def normalize_consensus(
    input_path: str,
    output_path: str,
    method: str = "median",
) -> int:
    """Normalize a consensus map across channels.

    Parameters
    ----------
    input_path : str
        Path to the input consensusXML file.
    output_path : str
        Path for the output normalized consensusXML file.
    method : str
        Normalization method: 'median' or 'quantile'.

    Returns
    -------
    int
        Number of consensus features in the normalized map.
    """
    if method not in VALID_METHODS:
        raise ValueError(f"Unknown method '{method}'. Valid: {VALID_METHODS}")

    cm = oms.ConsensusMap()
    oms.ConsensusXMLFile().load(input_path, cm)

    if method == "median":
        algo = oms.ConsensusMapNormalizerAlgorithmMedian()
        # method=0 corresponds to NM_SCALE (ratio-based normalization)
        algo.normalizeMaps(cm, 0, "", "")
    elif method == "quantile":
        algo = oms.ConsensusMapNormalizerAlgorithmQuantile()
        algo.normalizeMaps(cm)

    oms.ConsensusXMLFile().store(output_path, cm)
    return cm.size()


def create_synthetic_consensus_map(
    output_path: str,
    n_features: int = 5,
    n_channels: int = 2,
    intensity_bias: float = 2.0,
) -> None:
    """Create a synthetic consensusXML file with intensity bias in one channel.

    Parameters
    ----------
    output_path : str
        Path for the output consensusXML file.
    n_features : int
        Number of consensus features to generate.
    n_channels : int
        Number of channels (map indices).
    intensity_bias : float
        Multiplicative intensity bias applied to channel 1 relative to channel 0.
    """
    cm = oms.ConsensusMap()

    # Set up column headers
    headers = cm.getColumnHeaders()
    for i in range(n_channels):
        h = oms.ColumnHeader()
        h.filename = f"file{i}.mzML"
        headers[i] = h
    cm.setColumnHeaders(headers)

    for j in range(n_features):
        base_mz = 400.0 + j * 100.0
        base_rt = 50.0 + j * 30.0
        base_intensity = 10000.0 + j * 5000.0

        cf = oms.ConsensusFeature()
        cf.setMZ(base_mz)
        cf.setRT(base_rt)

        total_intensity = 0.0
        for ch in range(n_channels):
            p = oms.Peak2D()
            p.setMZ(base_mz)
            p.setRT(base_rt)
            if ch == 0:
                p.setIntensity(base_intensity)
            else:
                p.setIntensity(base_intensity * intensity_bias)
            total_intensity += p.getIntensity()
            cf.insert(ch, p, j)

        cf.setIntensity(total_intensity / n_channels)
        cm.push_back(cf)

    oms.ConsensusXMLFile().store(output_path, cm)


@click.command(help="Normalize intensity values in a consensusXML file.")
@click.option("--input", "input_path", required=True, help="Input consensusXML file")
@click.option("--output", "output_path", required=True, help="Output normalized consensusXML")
@click.option(
    "--method",
    default="median",
    type=click.Choice(VALID_METHODS),
    help="Normalization method (default: median)",
)
def main(input_path, output_path, method) -> None:
    n_features = normalize_consensus(input_path, output_path, method=method)
    click.echo(
        f"Normalized {n_features} features using '{method}' method, "
        f"saved to {output_path}"
    )


if __name__ == "__main__":
    main()
