"""
idXML to TSV Exporter
=====================
Export peptide and protein identifications from idXML to flat TSV format.

Usage
-----
    python idxml_to_tsv_exporter.py --input results.idXML --output results.tsv
"""

import csv
import sys
from typing import List

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def _make_peptide_id_list():
    """Create a peptide ID container compatible with the installed pyopenms version."""
    if hasattr(oms, "PeptideIdentificationList"):
        return oms.PeptideIdentificationList()
    return []


def load_idxml(input_path: str) -> tuple:
    """Load an idXML file and return (protein_ids, peptide_ids)."""
    protein_ids = []
    peptide_ids = _make_peptide_id_list()
    oms.IdXMLFile().load(input_path, protein_ids, peptide_ids)
    return protein_ids, list(peptide_ids)


def export_peptide_ids(peptide_ids: List[oms.PeptideIdentification], output_path: str) -> dict:
    """Export peptide identifications to TSV.

    Returns statistics about the export.
    """
    total_psms = 0

    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "spectrum_reference", "rt", "mz", "sequence", "charge",
            "score", "score_type", "protein_accessions",
        ])

        for pep_id in peptide_ids:
            spec_ref = pep_id.getMetaValue("spectrum_reference") if pep_id.metaValueExists(
                "spectrum_reference"
            ) else ""
            rt = pep_id.getRT()
            mz = pep_id.getMZ()
            score_type = pep_id.getScoreType()

            for hit in pep_id.getHits():
                accessions = ";".join(
                    ev.getProteinAccession() for ev in hit.getPeptideEvidences()
                )
                writer.writerow([
                    spec_ref,
                    f"{rt:.4f}",
                    f"{mz:.6f}",
                    hit.getSequence().toString(),
                    hit.getCharge(),
                    f"{hit.getScore():.6f}",
                    score_type,
                    accessions,
                ])
                total_psms += 1

    return {"peptide_ids": len(peptide_ids), "total_psms": total_psms}


def export_idxml(input_path: str, output_path: str) -> dict:
    """Export an idXML file to TSV format.

    Returns export statistics.
    """
    protein_ids, peptide_ids = load_idxml(input_path)
    stats = export_peptide_ids(peptide_ids, output_path)
    stats["protein_ids"] = len(protein_ids)
    return stats


def create_synthetic_idxml(output_path: str) -> None:
    """Create a synthetic idXML file for testing."""
    protein_id = oms.ProteinIdentification()
    protein_id.setSearchEngine("SEQUEST")
    protein_id.setScoreType("XCorr")
    protein_id.setIdentifier("run1")

    prot_hit = oms.ProteinHit()
    prot_hit.setAccession("P12345")
    prot_hit.setScore(100.0)
    protein_id.setHits([prot_hit])

    peptide_ids = _make_peptide_id_list()
    sequences = ["ACDEFGHIK", "MNPQRSTWY", "ACDEFGHIK"]
    for i, seq in enumerate(sequences):
        pep_id = oms.PeptideIdentification()
        pep_id.setRT(100.0 + i * 10)
        pep_id.setMZ(500.0 + i * 50)
        pep_id.setScoreType("XCorr")
        pep_id.setIdentifier("run1")

        pep_hit = oms.PeptideHit()
        pep_hit.setSequence(oms.AASequence.fromString(seq))
        pep_hit.setCharge(2)
        pep_hit.setScore(2.5 + i * 0.5)
        pep_hit.setRank(1)

        ev = oms.PeptideEvidence()
        ev.setProteinAccession("P12345")
        pep_hit.setPeptideEvidences([ev])

        pep_id.setHits([pep_hit])
        if hasattr(peptide_ids, "push_back"):
            peptide_ids.push_back(pep_id)
        else:
            peptide_ids.append(pep_id)

    oms.IdXMLFile().store(output_path, [protein_id], peptide_ids)


@click.command(help="Export idXML to flat TSV format.")
@click.option("--input", "input", required=True, help="Input idXML file")
@click.option("--output", required=True, help="Output TSV file")
def main(input, output) -> None:
    stats = export_idxml(input, output)
    print(f"Exported {stats['total_psms']} PSMs from {stats['peptide_ids']} spectra to {output}")


if __name__ == "__main__":
    main()
