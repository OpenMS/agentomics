"""
Posterior Error Probability Estimator
=====================================
Estimate posterior error probabilities (PEP) for peptide-spectrum matches
using a mixture model approach.

Usage
-----
    python posterior_error_probability.py --input peptides.idXML --output pep.idXML
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


def estimate_pep(input_path, output_path):
    """Estimate posterior error probabilities for PSMs.

    Extracts search engine scores, fits a two-component mixture model,
    and assigns PEP values to each PSM.

    Parameters
    ----------
    input_path : str
        Path to input idXML file.
    output_path : str
        Path to write output idXML with PEP scores.

    Returns
    -------
    int
        Number of PSMs scored.
    """
    protein_ids, peptide_ids = _load_idxml(input_path)

    # Extract all scores
    scores = []
    for pep_id in peptide_ids:
        for hit in pep_id.getHits():
            scores.append(hit.getScore())

    if not scores:
        oms.IdXMLFile().store(output_path, protein_ids, peptide_ids)
        return 0

    # Fit the PEP model
    pep_model = oms.PosteriorErrorProbabilityModel()
    pep_model.fit(scores, "GAUSS")

    # Apply PEP values back to hits (use index access since iteration yields copies)
    total_psms = 0
    for i in range(len(peptide_ids)):
        pep_id = peptide_ids[i]
        hits = pep_id.getHits()
        updated_hits = []
        for hit in hits:
            prob = pep_model.computeProbability(hit.getScore())
            hit.setScore(prob)
            updated_hits.append(hit)
            total_psms += 1
        pep_id.setHits(updated_hits)
        pep_id.setScoreType("Posterior Error Probability")
        pep_id.setHigherScoreBetter(False)
        peptide_ids[i] = pep_id

    oms.IdXMLFile().store(output_path, protein_ids, peptide_ids)
    return total_psms


@click.command(help="Estimate posterior error probabilities for PSMs.")
@click.option("--input", "input_path", required=True, help="Input idXML file")
@click.option("--output", "output_path", required=True, help="Output idXML with PEP scores")
def main(input_path, output_path):
    n = estimate_pep(input_path, output_path)
    print(f"PEP estimation complete. Scored {n} PSMs.")


if __name__ == "__main__":
    main()
