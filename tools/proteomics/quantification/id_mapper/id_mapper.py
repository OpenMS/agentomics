"""
ID Mapper
=========
Map peptide identifications onto features by matching RT and m/z values
within configurable tolerances.

Usage
-----
    python id_mapper.py --features features.featureXML --ids peptides.idXML \
        --output annotated.featureXML --rt-tol 5.0 --mz-tol 10.0
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def map_ids_to_features(
    features_path: str,
    ids_path: str,
    output_path: str,
    rt_tol: float = 5.0,
    mz_tol: float = 10.0,
) -> int:
    """Map peptide identifications to features based on RT and m/z proximity.

    Parameters
    ----------
    features_path : str
        Path to the input featureXML file.
    ids_path : str
        Path to the input idXML file.
    output_path : str
        Path for the output annotated featureXML file.
    rt_tol : float
        RT tolerance in seconds.
    mz_tol : float
        m/z tolerance in Da.

    Returns
    -------
    int
        Number of features with at least one peptide identification.
    """
    # Load feature map
    fm = oms.FeatureMap()
    oms.FeatureXMLFile().load(features_path, fm)

    # Load identifications
    protein_ids = []
    peptide_ids = oms.PeptideIdentificationList()
    oms.IdXMLFile().load(ids_path, protein_ids, peptide_ids)

    # Configure IDMapper
    mapper = oms.IDMapper()
    params = mapper.getDefaults()
    params.setValue("rt_tolerance", rt_tol)
    params.setValue("mz_tolerance", mz_tol)
    params.setValue("mz_measure", "Da")
    mapper.setParameters(params)

    # Annotate features with identifications
    # The FeatureMap overload requires: (fm, ids, prot_ids, use_centroid_rt, use_centroid_mz, spectra)
    empty_exp = oms.MSExperiment()
    mapper.annotate(fm, peptide_ids, protein_ids, True, True, empty_exp)

    # Count features with identifications
    n_annotated = 0
    for feature in fm:
        if feature.getPeptideIdentifications():
            n_annotated += 1

    # Store result
    oms.FeatureXMLFile().store(output_path, fm)

    return n_annotated


def create_synthetic_featurexml(
    output_path: str,
    mz: float = 500.0,
    rt: float = 100.0,
    intensity: float = 1e5,
) -> None:
    """Create a synthetic featureXML with a single feature.

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
    """
    fm = oms.FeatureMap()
    f = oms.Feature()
    f.setMZ(mz)
    f.setRT(rt)
    f.setIntensity(intensity)
    f.setCharge(2)
    f.setOverallQuality(0.9)
    fm.push_back(f)
    oms.FeatureXMLFile().store(output_path, fm)


def create_synthetic_idxml(
    output_path: str,
    mz: float = 500.0,
    rt: float = 100.0,
    sequence: str = "PEPTIDEK",
    score: float = 0.95,
) -> None:
    """Create a synthetic idXML with a single peptide identification.

    Parameters
    ----------
    output_path : str
        Path for the output idXML file.
    mz : float
        m/z of the peptide identification.
    rt : float
        RT of the peptide identification in seconds.
    sequence : str
        Peptide sequence string.
    score : float
        Score of the peptide hit.
    """
    protein_id = oms.ProteinIdentification()
    protein_id.setIdentifier("search_1")

    peptide_id = oms.PeptideIdentification()
    peptide_id.setRT(rt)
    peptide_id.setMZ(mz)
    peptide_id.setIdentifier("search_1")
    peptide_id.setScoreType("PEP")
    peptide_id.setHigherScoreBetter(True)

    hit = oms.PeptideHit()
    hit.setSequence(oms.AASequence.fromString(sequence))
    hit.setScore(score)
    hit.setCharge(2)
    peptide_id.setHits([hit])

    pid_list = oms.PeptideIdentificationList()
    pid_list.push_back(peptide_id)
    oms.IdXMLFile().store(output_path, [protein_id], pid_list)


@click.command(help="Map peptide identifications onto features.")
@click.option("--features", "features_path", required=True, help="Input featureXML file")
@click.option("--ids", "ids_path", required=True, help="Input idXML file")
@click.option("--output", "output_path", required=True, help="Output annotated featureXML")
@click.option("--rt-tol", default=5.0, help="RT tolerance in seconds (default: 5.0)")
@click.option("--mz-tol", default=10.0, help="m/z tolerance in Da (default: 10.0)")
def main(features_path, ids_path, output_path, rt_tol, mz_tol) -> None:
    n_annotated = map_ids_to_features(
        features_path, ids_path, output_path, rt_tol=rt_tol, mz_tol=mz_tol
    )
    click.echo(
        f"Mapped IDs to features: {n_annotated} features annotated, "
        f"saved to {output_path}"
    )


if __name__ == "__main__":
    main()
