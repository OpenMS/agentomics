"""
Protein In-Silico Digest
========================
Digest a protein sequence with a chosen enzyme and report the resulting
peptides together with their masses using pyopenms.

Supported enzymes (examples): Trypsin, Lys-C, Arg-C, Asp-N, Glu-C, Chymotrypsin
For the full list call:  python protein_digest.py --list-enzymes

Usage
-----
    python protein_digest.py --sequence MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELALGRFWDYLRWVQTLSEQVQEELLSSQVTQELRALMDETMKELKAYKSELEEQLTPVAEETRARLSKELQAAQARLGADVLASHGRLVQYRGEVQAMLGQSTEELRVRLASHLRKLRKRLLRDADDLQKRLAVYQAGAREGAERGLSAIRERLGPLVEQGRVRAATVGSLAGQPLQERAQAWGERLRARMEEMGSRTRDRLDEVKEQVAEVRAKLEEQAQQIRLQAEAFQARLKSWFEPLVEDMQRQWAGLVEKVQAAVGTSAAPVPSDNH --enzyme Trypsin
    python protein_digest.py --sequence MKVLWAALLVTFLAGC --enzyme Lys-C --missed-cleavages 2
"""

import argparse
import sys

try:
    import pyopenms as oms
except ImportError:
    sys.exit(
        "pyopenms is required. Install it with:  pip install pyopenms"
    )


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


def main():
    parser = argparse.ArgumentParser(
        description="In-silico protein digestion using pyopenms."
    )
    parser.add_argument(
        "--sequence",
        help="Single-letter amino acid sequence of the protein",
    )
    parser.add_argument(
        "--enzyme",
        default="Trypsin",
        help="Digestion enzyme name (default: Trypsin)",
    )
    parser.add_argument(
        "--missed-cleavages",
        type=int,
        default=0,
        dest="missed_cleavages",
        help="Maximum missed cleavages (default: 0)",
    )
    parser.add_argument(
        "--min-length",
        type=int,
        default=6,
        dest="min_length",
        help="Minimum peptide length (default: 6)",
    )
    parser.add_argument(
        "--max-length",
        type=int,
        default=40,
        dest="max_length",
        help="Maximum peptide length (default: 40)",
    )
    parser.add_argument(
        "--list-enzymes",
        action="store_true",
        dest="list_enzymes",
        help="List all available enzyme names and exit",
    )
    args = parser.parse_args()

    if args.list_enzymes:
        enzymes = list_enzymes()
        print("Available enzymes:")
        for name in enzymes:
            print(f"  {name}")
        return

    if not args.sequence:
        parser.error("--sequence is required unless --list-enzymes is used.")

    peptides = digest_protein(
        args.sequence,
        enzyme=args.enzyme,
        missed_cleavages=args.missed_cleavages,
        min_length=args.min_length,
        max_length=args.max_length,
    )

    print(
        f"Enzyme: {args.enzyme}  |  Missed cleavages ≤ {args.missed_cleavages}  "
        f"|  Length {args.min_length}–{args.max_length}"
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
