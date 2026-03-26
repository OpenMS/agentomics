"""
RT Alignment - Identification Based
====================================
Align retention times between runs using shared peptide identifications.

This tool uses peptide identifications present in both reference and input
to compute a retention time transformation, then applies it to correct
systematic RT shifts between LC-MS/MS runs.

Usage
-----
    python rt_alignment_identification.py \\
        --reference ref.idXML --input sample.idXML --output aligned.idXML --model linear
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

# Peptide sequences used for synthetic data generation
_SYNTHETIC_SEQUENCES = [
    "PEPTIDEK",
    "EDITPEPK",
    "ACDEKFGH",
    "KMLNPQRS",
    "TVWYKACDE",
    "AAELPKDR",
    "GFHKDTMR",
    "LNQYVDKW",
    "MPVSQETK",
    "RDFGHKLM",
]


def create_synthetic_idxml(
    output_path: str, n_ids: int = 5, rt_offset: float = 0.0
) -> None:
    """Create a synthetic idXML file for testing.

    Parameters
    ----------
    output_path : str
        Path to write the idXML file.
    n_ids : int
        Number of peptide identifications to generate (max 10).
    rt_offset : float
        Offset to add to retention times (seconds).
    """
    n_ids = min(n_ids, len(_SYNTHETIC_SEQUENCES))

    protein_ids = [oms.ProteinIdentification()]
    peptide_ids = oms.PeptideIdentificationList()

    for i in range(n_ids):
        pid = oms.PeptideIdentification()
        pid.setRT(100.0 + i * 100.0 + rt_offset)
        pid.setMZ(500.0 + i * 100.0)

        hit = oms.PeptideHit()
        hit.setSequence(oms.AASequence.fromString(_SYNTHETIC_SEQUENCES[i]))
        hit.setScore(0.95)
        pid.setHits([hit])

        peptide_ids.push_back(pid)

    oms.IdXMLFile().store(output_path, protein_ids, peptide_ids)


def _idxml_to_feature_map(
    protein_ids: list, peptide_ids: "oms.PeptideIdentificationList"
) -> "oms.FeatureMap":
    """Convert peptide identifications to a FeatureMap with embedded IDs.

    Parameters
    ----------
    protein_ids : list
        List of ProteinIdentification objects.
    peptide_ids : PeptideIdentificationList
        Peptide identifications to embed.

    Returns
    -------
    oms.FeatureMap
        FeatureMap with one feature per peptide identification.
    """
    fm = oms.FeatureMap()
    fm.setUniqueId(1)

    for i in range(peptide_ids.size()):
        pid = peptide_ids.at(i)
        f = oms.Feature()
        f.setMZ(pid.getMZ())
        f.setRT(pid.getRT())
        f.setIntensity(1e5)
        f.setCharge(2)
        f.setOverallQuality(1.0)
        f.setUniqueId(i + 2)

        pid_list = oms.PeptideIdentificationList()
        pid_list.push_back(pid)
        f.setPeptideIdentifications(pid_list)

        fm.push_back(f)

    return fm


def align_identification(
    reference_path: str,
    input_path: str,
    output_path: str,
    model: str = "linear",
) -> dict:
    """Align peptide identifications from input to reference using shared IDs.

    Parameters
    ----------
    reference_path : str
        Path to the reference idXML file.
    input_path : str
        Path to the input idXML file to be aligned.
    output_path : str
        Path to write the aligned idXML file.
    model : str
        Transformation model to fit: "linear" or "interpolated".

    Returns
    -------
    dict
        Statistics about the alignment:
        - n_peptide_ids: number of peptide IDs in the input
        - n_trafo_points: number of transformation data points
        - model: transformation model used
    """
    # Load reference idXML
    ref_prot = []
    ref_pep = oms.PeptideIdentificationList()
    oms.IdXMLFile().load(reference_path, ref_prot, ref_pep)

    # Load input idXML
    inp_prot = []
    inp_pep = oms.PeptideIdentificationList()
    oms.IdXMLFile().load(input_path, inp_prot, inp_pep)

    # Convert to FeatureMaps for the alignment algorithm
    ref_fm = _idxml_to_feature_map(ref_prot, ref_pep)
    inp_fm = _idxml_to_feature_map(inp_prot, inp_pep)

    # Set up and run identification-based alignment
    aligner = oms.MapAlignmentAlgorithmIdentification()
    params = aligner.getDefaults()
    params.setValue(b"min_run_occur", 2)
    aligner.setParameters(params)

    maps = [ref_fm, inp_fm]
    trafos = [oms.TransformationDescription(), oms.TransformationDescription()]

    aligner.setReference(ref_fm)
    aligner.align(maps, trafos, 0)

    # The transformation for the input is at index 1 (index 0 is reference)
    trafo = trafos[1]
    n_trafo_points = len(trafo.getDataPoints())

    # Fit a model to the transformation data points
    if n_trafo_points > 0:
        trafo.fitModel(model.encode())

    # Apply transformation to the input peptide identifications
    inp_pep_list = []
    for i in range(inp_pep.size()):
        inp_pep_list.append(inp_pep.at(i))

    aligned_pep = oms.PeptideIdentificationList()
    for pid in inp_pep_list:
        if n_trafo_points > 0:
            pid.setRT(trafo.apply(pid.getRT()))
        aligned_pep.push_back(pid)

    # Save aligned idXML
    out_prot = inp_prot if inp_prot else [oms.ProteinIdentification()]
    oms.IdXMLFile().store(output_path, out_prot, aligned_pep)

    return {
        "n_peptide_ids": inp_pep.size(),
        "n_trafo_points": n_trafo_points,
        "model": model,
    }


@click.command(help="Align retention times using shared peptide identifications.")
@click.option(
    "--reference", required=True, help="Reference idXML file."
)
@click.option(
    "--input", "input_path", required=True, help="Input idXML file to align."
)
@click.option(
    "--output", required=True, help="Output aligned idXML file."
)
@click.option(
    "--model",
    default="linear",
    type=click.Choice(["linear", "interpolated"]),
    help="Transformation model (default: linear).",
)
def main(reference: str, input_path: str, output: str, model: str) -> None:
    """CLI entry point for identification-based alignment."""
    stats = align_identification(reference, input_path, output, model=model)
    click.echo(
        f"Aligned {stats['n_peptide_ids']} peptide IDs "
        f"({stats['n_trafo_points']} transformation points, "
        f"model={stats['model']}) -> {output}"
    )


if __name__ == "__main__":
    main()
