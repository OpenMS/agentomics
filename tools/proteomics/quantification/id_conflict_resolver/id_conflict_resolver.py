"""
ID Conflict Resolver
====================
Resolve ambiguous peptide identification annotations on features,
keeping only the best-scoring hit per feature.

Usage
-----
    python id_conflict_resolver.py --input features.featureXML \
        --output resolved.featureXML
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def resolve_id_conflicts(input_path: str, output_path: str) -> int:
    """Resolve ambiguous peptide identifications on features.

    For each feature with multiple peptide identifications, only the
    best-scoring identification is retained. Removed identifications
    are moved to the unassigned peptide list.

    Parameters
    ----------
    input_path : str
        Path to the input featureXML file.
    output_path : str
        Path for the output resolved featureXML file.

    Returns
    -------
    int
        Number of features in the resolved map.
    """
    fm = oms.FeatureMap()
    oms.FeatureXMLFile().load(input_path, fm)

    resolver = oms.IDConflictResolverAlgorithm()
    resolver.resolve(fm)

    oms.FeatureXMLFile().store(output_path, fm)
    return fm.size()


def create_synthetic_conflicting_featurexml(
    output_path: str,
    mz: float = 500.0,
    rt: float = 100.0,
    intensity: float = 1e5,
    sequences: list = None,
    scores: list = None,
) -> None:
    """Create a featureXML with a feature having multiple conflicting peptide IDs.

    Parameters
    ----------
    output_path : str
        Path for the output featureXML file.
    mz : float
        m/z of the feature.
    rt : float
        RT of the feature in seconds.
    intensity : float
        Intensity of the feature.
    sequences : list, optional
        List of peptide sequences. Default: ["PEPTIDEK", "ANOTHERK"].
    scores : list, optional
        List of scores corresponding to each sequence. Default: [0.9, 0.5].
    """
    if sequences is None:
        sequences = ["PEPTIDEK", "ANOTHERK"]
    if scores is None:
        scores = [0.9, 0.5]

    fm = oms.FeatureMap()

    # Set up protein identification
    prot_id = oms.ProteinIdentification()
    prot_id.setIdentifier("search_1")
    prot_id.setScoreType("Mascot")
    prot_id.setHigherScoreBetter(True)
    fm.setProteinIdentifications([prot_id])

    f = oms.Feature()
    f.setMZ(mz)
    f.setRT(rt)
    f.setIntensity(intensity)
    f.setCharge(2)
    f.setOverallQuality(0.9)

    # Add multiple peptide identifications to the feature
    pid_list = oms.PeptideIdentificationList()
    for seq, score in zip(sequences, scores):
        pid = oms.PeptideIdentification()
        pid.setRT(rt)
        pid.setMZ(mz)
        pid.setIdentifier("search_1")
        pid.setScoreType("Mascot")
        pid.setHigherScoreBetter(True)

        hit = oms.PeptideHit()
        hit.setSequence(oms.AASequence.fromString(seq))
        hit.setScore(score)
        hit.setCharge(2)
        pid.setHits([hit])

        pid_list.push_back(pid)

    f.setPeptideIdentifications(pid_list)
    fm.push_back(f)

    oms.FeatureXMLFile().store(output_path, fm)


@click.command(help="Resolve ambiguous peptide identifications on features.")
@click.option("--input", "input_path", required=True, help="Input featureXML file")
@click.option("--output", "output_path", required=True, help="Output resolved featureXML")
def main(input_path, output_path) -> None:
    n_features = resolve_id_conflicts(input_path, output_path)
    click.echo(f"Resolved ID conflicts for {n_features} features, saved to {output_path}")


if __name__ == "__main__":
    main()
