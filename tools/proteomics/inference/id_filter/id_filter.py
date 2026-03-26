"""
ID Filter
=========
Filter peptide/protein identifications by score threshold and decoy status.

Usage
-----
    python id_filter.py --input peptides.idXML --output filtered.idXML --score-threshold 0.05
    python id_filter.py --input peptides.idXML --output filtered.idXML --remove-decoys
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


def _count_hits(peptide_ids):
    """Count total peptide hits across all identifications."""
    return sum(len(pid.getHits()) for pid in peptide_ids)


def filter_ids(input_path, output_path, score_threshold=0.05, remove_decoys=False):
    """Filter peptide identifications by score and/or decoy status.

    Parameters
    ----------
    input_path : str
        Path to input idXML file.
    output_path : str
        Path to write filtered idXML.
    score_threshold : float or None
        Score threshold for filtering (interpretation depends on score type).
    remove_decoys : bool
        If True, remove decoy hits.

    Returns
    -------
    dict
        Dictionary with 'before' and 'after' hit counts.
    """
    protein_ids, peptide_ids = _load_idxml(input_path)

    before_count = _count_hits(peptide_ids)

    id_filter = oms.IDFilter()

    if score_threshold is not None:
        id_filter.filterHitsByScore(peptide_ids, score_threshold)

    if remove_decoys:
        id_filter.removeDecoyHits(peptide_ids)

    # Remove empty identifications
    id_filter.removeEmptyIdentifications(peptide_ids)

    after_count = _count_hits(peptide_ids)

    oms.IdXMLFile().store(output_path, protein_ids, peptide_ids)

    return {
        "before": before_count,
        "after": after_count,
        "removed": before_count - after_count,
    }


@click.command(help="Filter peptide identifications by score and decoy status.")
@click.option("--input", "input_path", required=True, help="Input idXML file")
@click.option("--output", "output_path", required=True, help="Output filtered idXML")
@click.option(
    "--score-threshold", type=float, default=0.05,
    help="Score threshold for filtering (default: 0.05)",
)
@click.option(
    "--remove-decoys", is_flag=True, default=False,
    help="Remove decoy hits",
)
def main(input_path, output_path, score_threshold, remove_decoys):
    result = filter_ids(input_path, output_path, score_threshold, remove_decoys)
    print(
        f"Filtering complete. {result['before']} -> {result['after']} hits "
        f"({result['removed']} removed)"
    )


if __name__ == "__main__":
    main()
