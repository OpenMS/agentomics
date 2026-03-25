"""
FASTA In-Silico Digest Stats
=============================
Digest a FASTA database in silico and report peptide statistics including
peptide count, length distribution, mass distribution, and unique peptides.

Usage
-----
    python fasta_in_silico_digest_stats.py --input db.fasta --enzyme Trypsin --missed-cleavages 2 --output stats.tsv
"""

import csv
from typing import Dict, List

import click
import pyopenms as oms


def load_fasta(input_path: str) -> List[oms.FASTAEntry]:
    """Load entries from a FASTA file."""
    entries = []
    fasta_file = oms.FASTAFile()
    fasta_file.load(input_path, entries)
    return entries


def digest_fasta(
    input_path: str,
    enzyme: str = "Trypsin",
    missed_cleavages: int = 0,
    min_length: int = 6,
    max_length: int = 50,
) -> dict:
    """Digest all proteins in a FASTA and return peptide statistics.

    Returns a dict with:
    - protein_count: number of proteins digested
    - total_peptides: total peptide count
    - unique_peptides: unique peptide sequences
    - length_distribution: dict of length -> count
    - mass_stats: min, max, mean mass
    - peptides: list of dicts with sequence, mass, length, protein info
    """
    entries = load_fasta(input_path)
    digestor = oms.ProteaseDigestion()
    digestor.setEnzyme(enzyme)
    digestor.setMissedCleavages(missed_cleavages)

    all_peptides: List[dict] = []
    unique_seqs: set = set()

    for entry in entries:
        peptides: List[oms.AASequence] = []
        digestor.digest(oms.AASequence.fromString(entry.sequence), peptides)
        for pep in peptides:
            seq_str = pep.toString()
            seq_len = len(seq_str)
            if seq_len < min_length or seq_len > max_length:
                continue
            mass = pep.getMonoWeight()
            unique_seqs.add(seq_str)
            all_peptides.append({
                "sequence": seq_str,
                "length": seq_len,
                "mass": round(mass, 4),
                "protein": entry.identifier.split()[0],
            })

    # Length distribution
    length_dist: Dict[int, int] = {}
    masses: List[float] = []
    for pep in all_peptides:
        length_dist[pep["length"]] = length_dist.get(pep["length"], 0) + 1
        masses.append(pep["mass"])

    mass_stats = {}
    if masses:
        mass_stats = {
            "min": round(min(masses), 4),
            "max": round(max(masses), 4),
            "mean": round(sum(masses) / len(masses), 4),
        }

    return {
        "protein_count": len(entries),
        "enzyme": enzyme,
        "missed_cleavages": missed_cleavages,
        "total_peptides": len(all_peptides),
        "unique_peptides": len(unique_seqs),
        "length_distribution": dict(sorted(length_dist.items())),
        "mass_stats": mass_stats,
        "peptides": all_peptides,
    }


def write_tsv(stats: dict, output_path: str) -> None:
    """Write peptide digest results to a TSV file."""
    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["sequence", "length", "mass", "protein"])
        for pep in stats["peptides"]:
            writer.writerow([pep["sequence"], pep["length"], pep["mass"], pep["protein"]])


@click.command(help="Digest a FASTA database and report peptide statistics.")
@click.option("--input", "input", required=True, help="Input FASTA file")
@click.option("--enzyme", default="Trypsin", help="Enzyme name (default: Trypsin)")
@click.option("--missed-cleavages", type=int, default=0, help="Missed cleavages (default: 0)")
@click.option("--min-length", type=int, default=6, help="Min peptide length (default: 6)")
@click.option("--max-length", type=int, default=50, help="Max peptide length (default: 50)")
@click.option("--output", required=True, help="Output TSV file")
def main(input, enzyme, missed_cleavages, min_length, max_length, output) -> None:
    stats = digest_fasta(input, enzyme, missed_cleavages, min_length, max_length)
    write_tsv(stats, output)

    print(f"Proteins: {stats['protein_count']}")
    print(f"Total peptides: {stats['total_peptides']}")
    print(f"Unique peptides: {stats['unique_peptides']}")
    if stats["mass_stats"]:
        print(f"Mass range: {stats['mass_stats']['min']} - {stats['mass_stats']['max']}")
    print(f"Results written to {output}")


if __name__ == "__main__":
    main()
