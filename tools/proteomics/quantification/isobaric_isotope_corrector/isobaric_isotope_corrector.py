"""
Isobaric Isotope Corrector
===========================
Correct isotopic impurities in isobaric labeling quantification data.

Reads a consensusXML file produced by isobaric channel extraction and applies
the manufacturer-provided isotope correction matrix to remove cross-channel
contamination from isotopic impurities.

Usage
-----
    python isobaric_isotope_corrector.py --input quant.consensusXML --method tmt6plex --output corrected.consensusXML
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
    "itraq4plex": oms.ItraqFourPlexQuantitationMethod,
    "itraq8plex": oms.ItraqEightPlexQuantitationMethod,
}


def correct_isotope_impurities(
    input_path: str, method: str, output_path: str
) -> int:
    """Correct isotopic impurities in an isobaric quantification consensus map.

    Parameters
    ----------
    input_path : str
        Path to the input consensusXML file (from isobaric channel extraction).
    method : str
        Quantitation method name. One of: tmt6plex, tmt10plex,
        itraq4plex, itraq8plex.
    output_path : str
        Path for the output corrected consensusXML file.

    Returns
    -------
    int
        Number of consensus features in the corrected output.

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

    corrector = oms.IsobaricIsotopeCorrector()
    consensus_out = oms.ConsensusMap()
    corrector.correctIsotopicImpurities(consensus_in, consensus_out, quant_method)

    oms.ConsensusXMLFile().store(output_path, consensus_out)

    return consensus_out.size()


@click.command(
    help="Correct isotopic impurities in isobaric quantification data."
)
@click.option("--input", "input_path", required=True, help="Input consensusXML file.")
@click.option(
    "--method",
    required=True,
    type=click.Choice(sorted(SUPPORTED_METHODS), case_sensitive=False),
    help="Isobaric labeling method.",
)
@click.option(
    "--output", "output_path", required=True, help="Output corrected consensusXML file."
)
def main(input_path: str, method: str, output_path: str) -> None:
    n_features = correct_isotope_impurities(input_path, method, output_path)
    click.echo(f"Corrected {n_features} consensus features to {output_path}")


if __name__ == "__main__":
    main()
