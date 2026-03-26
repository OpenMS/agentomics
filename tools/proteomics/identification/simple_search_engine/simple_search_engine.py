"""
Simple Search Engine
====================
Run a simple peptide identification search against an mzML file and a
FASTA protein database using the OpenMS ``SimpleSearchEngineAlgorithm``.
Results are written in idXML format.

Usage
-----
    python simple_search_engine.py --input run.mzML --database proteins.fasta \\
        --output results.idXML --precursor-tol 10.0 --fragment-tol 20.0
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit(
        "pyopenms is required. Install it with:  pip install pyopenms"
    )


def search(
    input_path: str,
    database_path: str,
    output_path: str,
    precursor_tol: float = 10.0,
    fragment_tol: float = 20.0,
) -> int:
    """Run a simple peptide search engine.

    Parameters
    ----------
    input_path:
        Path to the input mzML file.
    database_path:
        Path to the FASTA protein database.
    output_path:
        Path to the output idXML file.
    precursor_tol:
        Precursor mass tolerance in ppm.
    fragment_tol:
        Fragment mass tolerance in ppm.

    Returns
    -------
    int
        Number of peptide-spectrum matches (PSMs) found.
    """
    sse = oms.SimpleSearchEngineAlgorithm()
    params = sse.getDefaults()
    params.setValue("precursor:mass_tolerance", float(precursor_tol))
    params.setValue("precursor:mass_tolerance_unit", "ppm")
    params.setValue("fragment:mass_tolerance", float(fragment_tol))
    params.setValue("fragment:mass_tolerance_unit", "ppm")
    sse.setParameters(params)

    protein_ids = []
    peptide_ids = oms.PeptideIdentificationList()
    sse.search(input_path, database_path, protein_ids, peptide_ids)

    # Count total PSMs (peptide hits across all identifications)
    psm_count = 0
    for pep_id in peptide_ids:
        psm_count += len(pep_id.getHits())

    # Write results to idXML
    oms.IdXMLFile().store(output_path, protein_ids, peptide_ids)

    return psm_count


@click.command(help="Run a simple peptide search engine on mzML data.")
@click.option("--input", "input_path", required=True, help="Input mzML file")
@click.option("--database", "database_path", required=True, help="FASTA protein database")
@click.option("--output", "output_path", required=True, help="Output idXML file")
@click.option(
    "--precursor-tol",
    type=float,
    default=10.0,
    help="Precursor mass tolerance in ppm (default: 10.0)",
)
@click.option(
    "--fragment-tol",
    type=float,
    default=20.0,
    help="Fragment mass tolerance in ppm (default: 20.0)",
)
def main(input_path, database_path, output_path, precursor_tol, fragment_tol):
    psm_count = search(
        input_path,
        database_path,
        output_path,
        precursor_tol=precursor_tol,
        fragment_tol=fragment_tol,
    )
    print(f"Found {psm_count} PSMs")
    print(f"Results written to {output_path}")


if __name__ == "__main__":
    main()
