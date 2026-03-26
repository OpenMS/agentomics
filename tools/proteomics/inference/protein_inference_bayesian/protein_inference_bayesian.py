"""
Bayesian Protein Inference
==========================
Infer proteins from peptide identifications using a Bayesian probabilistic
model from OpenMS.

Usage
-----
    python protein_inference_bayesian.py --input peptides.idXML --output proteins.idXML
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


def infer_proteins_bayesian(input_path, output_path):
    """Run Bayesian protein inference on peptide identifications.

    Parameters
    ----------
    input_path : str
        Path to input idXML file with peptide identifications.
    output_path : str
        Path to write output idXML with inferred proteins and posteriors.

    Returns
    -------
    int
        Number of proteins inferred.
    """
    protein_ids, peptide_ids = _load_idxml(input_path)

    algo = oms.BayesianProteinInferenceAlgorithm()
    algo.inferPosteriorProbabilities(protein_ids, peptide_ids, False)

    oms.IdXMLFile().store(output_path, protein_ids, peptide_ids)

    total_proteins = sum(len(pid.getHits()) for pid in protein_ids)
    return total_proteins


@click.command(help="Infer proteins using Bayesian probabilistic model.")
@click.option("--input", "input_path", required=True, help="Input idXML file")
@click.option("--output", "output_path", required=True, help="Output idXML with inferred proteins")
def main(input_path, output_path):
    n = infer_proteins_bayesian(input_path, output_path)
    print(f"Bayesian protein inference complete. {n} proteins inferred.")


if __name__ == "__main__":
    main()
