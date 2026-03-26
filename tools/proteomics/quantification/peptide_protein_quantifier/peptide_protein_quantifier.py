"""
Peptide and Protein Quantifier
==============================
Roll up peptide-level quantification to protein-level abundances from
an annotated consensusXML file using PeptideAndProteinQuant.

Usage
-----
    python peptide_protein_quantifier.py --input annotated.consensusXML \
        --output protein_quant.csv
"""

import csv
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def quantify_proteins(
    input_path: str,
    output_path: str,
    top_n: int = 3,
    include_all: bool = True,
) -> int:
    """Quantify proteins from an annotated consensus map.

    Parameters
    ----------
    input_path : str
        Path to the input consensusXML file (with peptide annotations).
    output_path : str
        Path for the output CSV file with protein abundances.
    top_n : int
        Number of top peptides to use for protein quantification.
    include_all : bool
        If True, include proteins with fewer than top_n peptides.

    Returns
    -------
    int
        Number of quantified proteins.
    """
    # Load consensus map
    cm = oms.ConsensusMap()
    oms.ConsensusXMLFile().load(input_path, cm)

    # Create experimental design from the consensus map
    ed = oms.ExperimentalDesign.fromConsensusMap(cm)

    # Set up quantifier
    quant = oms.PeptideAndProteinQuant()
    params = quant.getDefaults()
    params.setValue("top:N", top_n)
    params.setValue("top:include_all", "true" if include_all else "false")
    quant.setParameters(params)

    # Read quantitative data
    quant.readQuantData(cm, ed)

    # Quantify peptides
    empty_pid_list = oms.PeptideIdentificationList()
    quant.quantifyPeptides(empty_pid_list)

    # Quantify proteins
    prot_ids = cm.getProteinIdentifications()
    if prot_ids:
        quant.quantifyProteins(prot_ids[0])

    # Get statistics
    stats = quant.getStatistics()

    # Write results to CSV
    _write_results_csv(output_path, prot_ids, stats)

    return stats.quant_proteins


def _write_results_csv(
    output_path: str,
    prot_ids: list,
    stats,
) -> None:
    """Write quantification results to a CSV file.

    Parameters
    ----------
    output_path : str
        Path for the output CSV file.
    prot_ids : list
        List of ProteinIdentification objects.
    stats : PeptideAndProteinQuant_Statistics
        Quantification statistics.
    """
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "protein_accession",
            "score",
        ])

        if prot_ids:
            for hit in prot_ids[0].getHits():
                writer.writerow([
                    hit.getAccession(),
                    f"{hit.getScore():.6f}",
                ])

        # Write summary row
        writer.writerow([])
        writer.writerow(["# Summary"])
        writer.writerow(["total_features", stats.total_features])
        writer.writerow(["quant_peptides", stats.quant_peptides])
        writer.writerow(["quant_proteins", stats.quant_proteins])


def create_synthetic_annotated_consensus(
    output_path: str,
    proteins: dict = None,
) -> None:
    """Create a synthetic consensusXML with peptide annotations from multiple proteins.

    Parameters
    ----------
    output_path : str
        Path for the output consensusXML file.
    proteins : dict, optional
        Mapping of protein accession to list of (sequence, mz, rt, intensity) tuples.
        If None, creates a default 2-protein dataset.
    """
    if proteins is None:
        proteins = {
            "ProteinA": [
                ("PEPTIDEK", 500.0, 100.0, 10000.0),
                ("ANOTHERSEQK", 600.0, 200.0, 20000.0),
            ],
            "ProteinB": [
                ("DIFFERENTK", 700.0, 300.0, 15000.0),
                ("YETMORER", 800.0, 400.0, 25000.0),
            ],
        }

    cm = oms.ConsensusMap()

    # Set up column headers
    headers = cm.getColumnHeaders()
    h0 = oms.ColumnHeader()
    h0.filename = "file0.mzML"
    h0.label = "label0"
    headers[0] = h0
    cm.setColumnHeaders(headers)

    feature_idx = 0
    for prot_acc, peptides in proteins.items():
        for seq, mz, rt, intensity in peptides:
            cf = oms.ConsensusFeature()
            cf.setMZ(mz)
            cf.setRT(rt)
            cf.setIntensity(intensity)

            p = oms.Peak2D()
            p.setMZ(mz)
            p.setRT(rt)
            p.setIntensity(intensity)
            cf.insert(0, p, feature_idx)

            # Add peptide identification
            pid = oms.PeptideIdentification()
            pid.setRT(rt)
            pid.setMZ(mz)
            pid.setIdentifier("search_1")
            pid.setScoreType("PEP")
            pid.setHigherScoreBetter(False)

            hit = oms.PeptideHit()
            hit.setSequence(oms.AASequence.fromString(seq))
            hit.setScore(0.01)
            hit.setCharge(2)

            ev = oms.PeptideEvidence()
            ev.setProteinAccession(prot_acc)
            hit.setPeptideEvidences([ev])

            pid.setHits([hit])
            pid_list = oms.PeptideIdentificationList()
            pid_list.push_back(pid)
            cf.setPeptideIdentifications(pid_list)

            cm.push_back(cf)
            feature_idx += 1

    # Set protein identifications on the map
    prot_id = oms.ProteinIdentification()
    prot_id.setIdentifier("search_1")
    prot_id.setScoreType("PEP")
    prot_id.setHigherScoreBetter(False)
    for acc in proteins:
        ph = oms.ProteinHit()
        ph.setAccession(acc)
        ph.setScore(0.01)
        prot_id.insertHit(ph)
    cm.setProteinIdentifications([prot_id])

    oms.ConsensusXMLFile().store(output_path, cm)


@click.command(help="Quantify proteins from an annotated consensusXML file.")
@click.option("--input", "input_path", required=True, help="Input annotated consensusXML")
@click.option("--output", "output_path", required=True, help="Output protein quantification CSV")
@click.option("--top-n", default=3, help="Top N peptides for protein quant (default: 3)")
@click.option(
    "--include-all/--no-include-all",
    default=True,
    help="Include proteins with < top_n peptides (default: True)",
)
def main(input_path, output_path, top_n, include_all) -> None:
    n_proteins = quantify_proteins(
        input_path, output_path, top_n=top_n, include_all=include_all
    )
    click.echo(f"Quantified {n_proteins} proteins, saved to {output_path}")


if __name__ == "__main__":
    main()
