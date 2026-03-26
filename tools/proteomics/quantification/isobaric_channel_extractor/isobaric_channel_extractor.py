"""
Isobaric Channel Extractor
===========================
Extract reporter-ion channels from isobaric labeling experiments (TMT, iTRAQ).

Reads an mzML file containing MS2 spectra with reporter-ion peaks and extracts
per-channel intensities into a ConsensusMap stored as consensusXML.

Usage
-----
    python isobaric_channel_extractor.py --input run.mzML --method tmt6plex --output quant.consensusXML
    python isobaric_channel_extractor.py --input run.mzML --method itraq4plex --output quant.consensusXML
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


def extract_channels(input_path: str, method: str, output_path: str) -> int:
    """Extract isobaric reporter-ion channels from an mzML file.

    Parameters
    ----------
    input_path : str
        Path to the input mzML file containing MS2 spectra with reporter ions.
    method : str
        Quantitation method name. One of: tmt6plex, tmt10plex, tmt16plex,
        itraq4plex, itraq8plex.
    output_path : str
        Path for the output consensusXML file.

    Returns
    -------
    int
        Number of consensus features (MS2 scans) extracted.

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

    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    extractor = oms.IsobaricChannelExtractor(quant_method)
    consensus = oms.ConsensusMap()
    extractor.extractChannels(exp, consensus)

    oms.ConsensusXMLFile().store(output_path, consensus)

    return consensus.size()


@click.command(help="Extract isobaric reporter-ion channels from an mzML file.")
@click.option("--input", "input_path", required=True, help="Input mzML file.")
@click.option(
    "--method",
    required=True,
    type=click.Choice(sorted(SUPPORTED_METHODS), case_sensitive=False),
    help="Isobaric labeling method.",
)
@click.option("--output", "output_path", required=True, help="Output consensusXML file.")
def main(input_path: str, method: str, output_path: str) -> None:
    n_features = extract_channels(input_path, method, output_path)
    click.echo(f"Extracted {n_features} consensus features to {output_path}")


if __name__ == "__main__":
    main()
