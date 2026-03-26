"""
Isobaric Quantifier
====================
Quantify isobaric labeling experiments (TMT, iTRAQ) from a consensus map.

Reads a consensusXML file (from isobaric channel extraction or isotope
correction) and performs quantification with optional isotope correction
and normalization.

Usage
-----
    python isobaric_quantifier.py --input quant.consensusXML --method tmt6plex --output quantified.consensusXML
    python isobaric_quantifier.py --input quant.consensusXML --method itraq4plex --output quantified.consensusXML
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

SUPPORTED_METHODS = {
    "tmt6plex": oms.TMTSixPlexQuantitationMethod,
    "tmt10plex": oms.TMTTenPlexQuantitationMethod,
    "tmt16plex": oms.TMTSixteenPlexQuantitationMethod,
    "itraq4plex": oms.ItraqFourPlexQuantitationMethod,
    "itraq8plex": oms.ItraqEightPlexQuantitationMethod,
}


def quantify_isobaric(input_path: str, method: str, output_path: str) -> int:
    """Quantify isobaric labeling channels from a consensus map.

    Parameters
    ----------
    input_path : str
        Path to the input consensusXML file (from channel extraction).
    method : str
        Quantitation method name. One of: tmt6plex, tmt10plex, tmt16plex,
        itraq4plex, itraq8plex.
    output_path : str
        Path for the output quantified consensusXML file.

    Returns
    -------
    int
        Number of consensus features in the quantified output.

    Raises
    ------
    ValueError
        If the method name is not supported.
    """
    if method not in SUPPORTED_METHODS:
        raise ValueError(
            f"Unsupported method '{method}'. "
            f"Choose from: {', '.join(sorted(SUPPORTED_METHODS))}"
        )

    quant_method = SUPPORTED_METHODS[method]()

    consensus_in = oms.ConsensusMap()
    oms.ConsensusXMLFile().load(input_path, consensus_in)

    quantifier = oms.IsobaricQuantifier(quant_method)
    consensus_out = oms.ConsensusMap()
    quantifier.quantify(consensus_in, consensus_out)

    oms.ConsensusXMLFile().store(output_path, consensus_out)

    return consensus_out.size()


@click.command(help="Quantify isobaric labeling experiments from a consensus map.")
@click.option("--input", "input_path", required=True, help="Input consensusXML file.")
@click.option(
    "--method",
    required=True,
    type=click.Choice(sorted(SUPPORTED_METHODS), case_sensitive=False),
    help="Isobaric labeling method.",
)
@click.option(
    "--output", "output_path", required=True, help="Output quantified consensusXML file."
)
def main(input_path: str, method: str, output_path: str) -> None:
    n_features = quantify_isobaric(input_path, method, output_path)
    click.echo(f"Quantified {n_features} consensus features to {output_path}")


if __name__ == "__main__":
    main()
