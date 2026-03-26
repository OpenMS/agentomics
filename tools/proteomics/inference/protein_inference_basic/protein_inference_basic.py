"""
Basic Protein Inference
=======================
Infer proteins from peptide identifications using the basic aggregation
algorithm from OpenMS.

Usage
-----
    python protein_inference_basic.py --input peptides.idXML --output proteins.idXML
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


def infer_proteins(input_path, output_path):
    """Run basic protein inference on peptide identifications.

    Parameters
    ----------
    input_path : str
        Path to input idXML file with peptide identifications.
    output_path : str
        Path to write output idXML with inferred proteins.

    Returns
    -------
    int
        Number of proteins inferred.
    """
    protein_ids, peptide_ids = _load_idxml(input_path)

    algo = oms.BasicProteinInferenceAlgorithm()
    algo.run(peptide_ids, protein_ids)

    oms.IdXMLFile().store(output_path, protein_ids, peptide_ids)

    total_proteins = sum(len(pid.getHits()) for pid in protein_ids)
    return total_proteins


@click.command(help="Infer proteins from peptide identifications.")
@click.option("--input", "input_path", required=True, help="Input idXML file")
@click.option("--output", "output_path", required=True, help="Output idXML with inferred proteins")
def main(input_path, output_path):
    n = infer_proteins(input_path, output_path)
    print(f"Protein inference complete. {n} proteins inferred.")


if __name__ == "__main__":
    main()
