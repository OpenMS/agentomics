"""
Protein In-Silico Digest
========================
Digest a protein sequence with a chosen enzyme and report the resulting
peptides together with their masses using pyopenms.

Supported enzymes (examples): Trypsin, Lys-C, Arg-C, Asp-N, Glu-C, Chymotrypsin
For the full list call:  python protein_digest.py --list-enzymes

Usage
-----
    python protein_digest.py --sequence MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELAL \
        --enzyme Trypsin
    python protein_digest.py --sequence MKVLWAALLVTFLAGC --enzyme Lys-C --missed-cleavages 2
"""


import click
import pyopenms as oms


def list_enzymes() -> list:
    """Return all enzyme names registered in the pyopenms ProteaseDB."""
    db = oms.ProteaseDB()
    raw_names = []
    db.getAllNames(raw_names)
    return sorted(
        n.decode() if isinstance(n, bytes) else n for n in raw_names
    )


def digest_protein(
    sequence: str,
    enzyme: str = "Trypsin",
    missed_cleavages: int = 0,
    min_length: int = 6,
    max_length: int = 40,
) -> list:
    """Digest a protein sequence in silico.

    Parameters
    ----------
    sequence:
        Single-letter amino acid sequence of the protein.
    enzyme:
        Enzyme name as known to pyopenms ProteaseDB (default: ``"Trypsin"``).
    missed_cleavages:
        Maximum number of missed cleavages allowed (default: 0).
    min_length:
        Minimum peptide length to include (default: 6).
    max_length:
        Maximum peptide length to include (default: 40).

    Returns
    -------
    list of dict
        Each entry contains ``sequence``, ``start``, ``end``,
        ``monoisotopic_mass``, and ``missed_cleavages``.
    """
    protein_seq = oms.AASequence.fromString(sequence)
    digest = oms.ProteaseDigestion()
    digest.setEnzyme(enzyme)
    digest.setMissedCleavages(missed_cleavages)

    peptides = []
    digest.digest(protein_seq, peptides, min_length, max_length)

    results = []
    for pep in peptides:
        pep_str = pep.toString()
        mass = pep.getMonoWeight()
        results.append(
            {
                "sequence": pep_str,
                "length": pep.size(),
                "monoisotopic_mass": mass,
            }
        )
    return results


@click.command(help="In-silico protein digestion using pyopenms.")
@click.option("--sequence", default=None, help="Single-letter amino acid sequence of the protein")
@click.option("--enzyme", default="Trypsin", help="Digestion enzyme name (default: Trypsin)")
@click.option("--missed-cleavages", type=int, default=0, help="Maximum missed cleavages (default: 0)")
@click.option("--min-length", type=int, default=6, help="Minimum peptide length (default: 6)")
@click.option("--max-length", type=int, default=40, help="Maximum peptide length (default: 40)")
@click.option("--list-enzymes", "show_enzymes", is_flag=True, help="List all available enzyme names and exit")
def main(sequence, enzyme, missed_cleavages, min_length, max_length, show_enzymes):
    if show_enzymes:
        enzymes_list = list_enzymes()
        print("Available enzymes:")
        for name in enzymes_list:
            print(f"  {name}")
        return

    if not sequence:
        raise click.UsageError("--sequence is required unless --list-enzymes is used.")

    peptides = digest_protein(
        sequence,
        enzyme=enzyme,
        missed_cleavages=missed_cleavages,
        min_length=min_length,
        max_length=max_length,
    )

    print(
        f"Enzyme: {enzyme}  |  Missed cleavages \u2264 {missed_cleavages}  "
        f"|  Length {min_length}\u2013{max_length}"
    )
    print(f"Total peptides: {len(peptides)}\n")
    print(f"{'#':>4}  {'Sequence':<40}  {'Length':>6}  {'Mono Mass (Da)':>14}")
    print("-" * 72)
    for i, pep in enumerate(peptides, 1):
        print(
            f"{i:>4}  {pep['sequence']:<40}  {pep['length']:>6}  "
            f"{pep['monoisotopic_mass']:>14.6f}"
        )


if __name__ == "__main__":
    main()
