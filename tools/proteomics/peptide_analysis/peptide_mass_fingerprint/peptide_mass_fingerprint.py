"""
Peptide Mass Fingerprint
========================
Generate and match peptide mass fingerprints from FASTA databases.
Digests a protein with a specified enzyme and reports theoretical
peptide masses for fingerprint matching.

Usage
-----
    python peptide_mass_fingerprint.py --fasta db.fasta --accession P12345 --enzyme Trypsin --output fingerprint.tsv
"""

import csv
import sys
from typing import List

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit(
        "pyopenms is required. Install it with:  pip install pyopenms"
    )

PROTON = 1.007276


def load_fasta(fasta_path: str) -> dict:
    """Load FASTA file and return accession to sequence mapping.

    Parameters
    ----------
    fasta_path:
        Path to FASTA file.

    Returns
    -------
    dict
        Mapping of accession to protein sequence.
    """
    entries = []
    fasta_file = oms.FASTAFile()
    fasta_file.load(fasta_path, entries)
    proteins = {}
    for entry in entries:
        proteins[entry.identifier] = entry.sequence
    return proteins


def generate_fingerprint(
    protein_sequence: str,
    enzyme: str = "Trypsin",
    missed_cleavages: int = 1,
    min_mass: float = 500.0,
    max_mass: float = 4000.0,
) -> List[dict]:
    """Generate peptide mass fingerprint for a protein.

    Parameters
    ----------
    protein_sequence:
        Full protein amino acid sequence.
    enzyme:
        Enzyme name for digestion.
    missed_cleavages:
        Number of allowed missed cleavages.
    min_mass:
        Minimum peptide mass to include (Da).
    max_mass:
        Maximum peptide mass to include (Da).

    Returns
    -------
    list
        List of dicts with peptide sequence, mass, and position info.
    """
    digestion = oms.ProteaseDigestion()
    digestion.setEnzyme(enzyme)
    digestion.setMissedCleavages(missed_cleavages)

    aa_protein = oms.AASequence.fromString(protein_sequence)
    peptides = []
    digestion.digest(aa_protein, peptides)

    results = []
    for pep in peptides:
        pep_str = str(pep)
        mono_mass = pep.getMonoWeight()
        if mono_mass < min_mass or mono_mass > max_mass:
            continue

        # Find position in protein
        pos = protein_sequence.find(pep_str)

        results.append({
            "sequence": pep_str,
            "monoisotopic_mass": round(mono_mass, 6),
            "mz_1": round((mono_mass + PROTON) / 1, 6),
            "mz_2": round((mono_mass + 2 * PROTON) / 2, 6),
            "length": len(pep_str),
            "position": pos if pos >= 0 else "N/A",
        })

    results.sort(key=lambda x: x["monoisotopic_mass"])
    return results


def match_fingerprint(
    fingerprint: List[dict],
    observed_masses: List[float],
    tolerance_ppm: float = 10.0,
) -> List[dict]:
    """Match observed masses against a theoretical fingerprint.

    Parameters
    ----------
    fingerprint:
        Theoretical fingerprint from generate_fingerprint().
    observed_masses:
        List of observed monoisotopic masses.
    tolerance_ppm:
        Mass tolerance in ppm.

    Returns
    -------
    list
        List of matched peptide dicts with observed mass and ppm error.
    """
    matches = []
    for obs in observed_masses:
        for pep in fingerprint:
            theo = pep["monoisotopic_mass"]
            ppm_error = abs(obs - theo) / theo * 1e6
            if ppm_error <= tolerance_ppm:
                match = dict(pep)
                match["observed_mass"] = round(obs, 6)
                match["ppm_error"] = round(ppm_error, 2)
                matches.append(match)
    return matches


def write_tsv(records: List[dict], output_path: str) -> None:
    """Write fingerprint records to TSV.

    Parameters
    ----------
    records:
        List of record dicts.
    output_path:
        Output file path.
    """
    if not records:
        return
    fieldnames = list(records[0].keys())
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in records:
            writer.writerow(row)


@click.command(help="Generate/match peptide mass fingerprints from FASTA.")
@click.option("--fasta", required=True, help="FASTA database file")
@click.option("--accession", required=True, help="Protein accession to fingerprint")
@click.option("--enzyme", default="Trypsin", help="Enzyme name (default: Trypsin)")
@click.option("--missed-cleavages", type=int, default=1, help="Missed cleavages (default: 1)")
@click.option("--min-mass", type=float, default=500.0, help="Min peptide mass (default: 500)")
@click.option("--max-mass", type=float, default=4000.0, help="Max peptide mass (default: 4000)")
@click.option("--output", default=None, help="Output TSV file path")
def main(fasta, accession, enzyme, missed_cleavages, min_mass, max_mass, output):
    proteins = load_fasta(fasta)
    if accession not in proteins:
        sys.exit(f"Accession '{accession}' not found in FASTA file.")

    protein_seq = proteins[accession]
    print(f"Protein {accession}: {len(protein_seq)} amino acids")

    fingerprint = generate_fingerprint(
        protein_seq, enzyme=enzyme,
        missed_cleavages=missed_cleavages,
        min_mass=min_mass, max_mass=max_mass,
    )
    print(f"Generated {len(fingerprint)} peptide masses")

    for pep in fingerprint[:10]:
        print(f"  {pep['monoisotopic_mass']:.4f}  {pep['sequence']}")
    if len(fingerprint) > 10:
        print(f"  ... and {len(fingerprint) - 10} more")

    if output:
        write_tsv(fingerprint, output)
        print(f"\nFingerprint written to {output}")


if __name__ == "__main__":
    main()
