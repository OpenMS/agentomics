"""
Consensus ID
=============
Merge peptide identifications from multiple search engines using consensus
scoring algorithms.

Usage
-----
    python consensus_id.py --inputs search1.idXML search2.idXML --output consensus.idXML --algorithm best
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


ALGORITHMS = {
    "best": "ConsensusIDAlgorithmBest",
    "average": "ConsensusIDAlgorithmAverage",
    "ranks": "ConsensusIDAlgorithmRanks",
}


def _load_idxml(input_path):
    """Load an idXML file and return (protein_ids, peptide_ids)."""
    protein_ids = []
    peptide_ids = oms.PeptideIdentificationList()
    oms.IdXMLFile().load(input_path, protein_ids, peptide_ids)
    return protein_ids, peptide_ids


def consensus_id(input_paths, output_path, algorithm="best"):
    """Merge peptide identifications from multiple search engines.

    Parameters
    ----------
    input_paths : list of str
        Paths to input idXML files.
    output_path : str
        Path to write consensus idXML.
    algorithm : str
        Consensus algorithm: 'best', 'average', or 'ranks'.

    Returns
    -------
    int
        Number of consensus peptide identifications.
    """
    if algorithm not in ALGORITHMS:
        raise ValueError(
            f"Unknown algorithm '{algorithm}'. Choose from: {list(ALGORITHMS.keys())}"
        )

    # Load all input files and merge peptide IDs
    all_protein_ids = []
    all_peptide_ids = oms.PeptideIdentificationList()
    num_runs = 0

    for path in input_paths:
        protein_ids, peptide_ids = _load_idxml(path)
        if not all_protein_ids:
            all_protein_ids = protein_ids
        for pep_id in peptide_ids:
            all_peptide_ids.push_back(pep_id)
        num_runs += 1

    # Get the algorithm class
    algo_class = getattr(oms, ALGORITHMS[algorithm])
    algo = algo_class()

    # Apply consensus algorithm (requires number_of_runs)
    algo.apply(all_peptide_ids, num_runs)

    # Remove empty identifications after consensus
    id_filter = oms.IDFilter()
    id_filter.removeEmptyIdentifications(all_peptide_ids)

    # Fix identifiers: consensus may produce empty identifiers on merged IDs
    if all_protein_ids:
        run_id = all_protein_ids[0].getIdentifier()
        for i in range(len(all_peptide_ids)):
            pep_id = all_peptide_ids[i]
            if not pep_id.getIdentifier():
                pep_id.setIdentifier(run_id)
                all_peptide_ids[i] = pep_id

    oms.IdXMLFile().store(output_path, all_protein_ids, all_peptide_ids)

    return len(all_peptide_ids)


@click.command(help="Merge peptide identifications using consensus scoring.")
@click.option(
    "--inputs", "input_paths", required=True, multiple=True,
    help="Input idXML files (specify multiple times)",
)
@click.option("--output", "output_path", required=True, help="Output consensus idXML")
@click.option(
    "--algorithm", type=click.Choice(["best", "average", "ranks"]),
    default="best", help="Consensus algorithm (default: best)",
)
def main(input_paths, output_path, algorithm):
    n = consensus_id(list(input_paths), output_path, algorithm)
    print(f"Consensus ID complete. {n} peptide identifications in output.")


if __name__ == "__main__":
    main()
