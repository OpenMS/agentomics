"""
Peptide Indexer
===============
Map peptide identifications to proteins in a FASTA database using
PeptideIndexing from OpenMS.

Usage
-----
    python peptide_indexer.py --ids peptides.idXML --fasta database.fasta --output indexed.idXML
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def _load_idxml(input_path):
    """Load an idXML file and return (protein_ids, peptide_ids)."""
    protein_ids = []
    peptide_ids = oms.PeptideIdentificationList()
    oms.IdXMLFile().load(input_path, protein_ids, peptide_ids)
    return protein_ids, peptide_ids


def index_peptides(ids_path, fasta_path, output_path):
    """Map peptide identifications to proteins in a FASTA database.

    Parameters
    ----------
    ids_path : str
        Path to input idXML file with peptide identifications.
    fasta_path : str
        Path to FASTA database file.
    output_path : str
        Path to write output idXML with protein references.

    Returns
    -------
    int
        Number of peptide hits indexed.
    """
    protein_ids, peptide_ids = _load_idxml(ids_path)

    # Load FASTA entries
    fasta_entries = []
    oms.FASTAFile().load(fasta_path, fasta_entries)

    indexer = oms.PeptideIndexing()

    # Set parameters for more permissive matching
    params = indexer.getParameters()
    params.setValue("missing_decoy_action", "warn")
    params.setValue("enzyme:name", "no cleavage")
    params.setValue("enzyme:specificity", "none")
    indexer.setParameters(params)

    indexer.run(fasta_entries, protein_ids, peptide_ids)

    oms.IdXMLFile().store(output_path, protein_ids, peptide_ids)

    total_hits = sum(len(pid.getHits()) for pid in peptide_ids)
    return total_hits


@click.command(help="Map peptide identifications to proteins in a FASTA database.")
@click.option("--ids", "ids_path", required=True, help="Input idXML file")
@click.option("--fasta", "fasta_path", required=True, help="FASTA database file")
@click.option("--output", "output_path", required=True, help="Output idXML with protein references")
def main(ids_path, fasta_path, output_path):
    n = index_peptides(ids_path, fasta_path, output_path)
    print(f"Peptide indexing complete. {n} peptide hits indexed.")


if __name__ == "__main__":
    main()
