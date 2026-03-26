"""
FDR Estimator
=============
Estimate false discovery rates for peptide and protein identifications
using the target-decoy approach.

Usage
-----
    python fdr_estimator.py --input peptides.idXML --output fdr.idXML
    python fdr_estimator.py --input peptides.idXML --output fdr.idXML --protein
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


def _count_at_thresholds(peptide_ids, thresholds=(0.01, 0.05, 0.1)):
    """Count PSMs passing various FDR thresholds based on q-value score."""
    counts = {}
    for thresh in thresholds:
        n = 0
        for pep_id in peptide_ids:
            for hit in pep_id.getHits():
                score = hit.getScore()
                if score <= thresh:
                    n += 1
        counts[thresh] = n
    return counts


def estimate_fdr(input_path, output_path, protein_level=False):
    """Estimate FDR on peptide (and optionally protein) identifications.

    Parameters
    ----------
    input_path : str
        Path to input idXML file.
    output_path : str
        Path to write output idXML with q-values.
    protein_level : bool
        If True, also compute protein-level FDR.

    Returns
    -------
    dict
        Counts of PSMs at various FDR thresholds.
    """
    protein_ids, peptide_ids = _load_idxml(input_path)

    fdr = oms.FalseDiscoveryRate()
    fdr.apply(peptide_ids)

    if protein_level and protein_ids:
        fdr.apply(protein_ids)

    oms.IdXMLFile().store(output_path, protein_ids, peptide_ids)

    result = {
        "total_psms": sum(len(pid.getHits()) for pid in peptide_ids),
        "counts_at_fdr": _count_at_thresholds(peptide_ids),
    }
    return result


@click.command(help="Estimate FDR for peptide/protein identifications.")
@click.option("--input", "input_path", required=True, help="Input idXML file")
@click.option("--output", "output_path", required=True, help="Output idXML with q-values")
@click.option("--protein/--no-protein", default=False, help="Also compute protein-level FDR")
def main(input_path, output_path, protein):
    result = estimate_fdr(input_path, output_path, protein_level=protein)
    total = result["total_psms"]
    counts = result["counts_at_fdr"]
    print(f"FDR estimation complete. Total PSMs: {total}")
    for thresh, n in sorted(counts.items()):
        print(f"  PSMs at {thresh:.0%} FDR: {n}")


if __name__ == "__main__":
    main()
