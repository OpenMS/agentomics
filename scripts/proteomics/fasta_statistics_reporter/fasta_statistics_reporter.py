"""
FASTA Statistics Reporter
=========================
Report statistics from a FASTA database: protein count, sequence lengths,
amino acid frequencies, tryptic peptide counts, and duplicate detection.

Usage
-----
    python fasta_statistics_reporter.py --input db.fasta --enzyme Trypsin --output stats.json
"""

import argparse
import json
import sys
from collections import Counter
from typing import Dict, List, Optional

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def load_fasta(input_path: str) -> List[oms.FASTAEntry]:
    """Load entries from a FASTA file."""
    entries = []
    fasta_file = oms.FASTAFile()
    fasta_file.load(input_path, entries)
    return entries


def compute_length_stats(entries: List[oms.FASTAEntry]) -> dict:
    """Compute min, max, mean, median sequence lengths."""
    lengths = sorted(len(e.sequence) for e in entries)
    if not lengths:
        return {"min": 0, "max": 0, "mean": 0.0, "median": 0.0}
    n = len(lengths)
    median = (lengths[n // 2] + lengths[(n - 1) // 2]) / 2.0
    return {
        "min": lengths[0],
        "max": lengths[-1],
        "mean": round(sum(lengths) / n, 2),
        "median": median,
    }


def compute_aa_frequency(entries: List[oms.FASTAEntry]) -> Dict[str, int]:
    """Count amino acid frequencies across all sequences."""
    counter: Counter = Counter()
    for entry in entries:
        counter.update(entry.sequence)
    return dict(sorted(counter.items()))


def count_tryptic_peptides(entries: List[oms.FASTAEntry], enzyme: str, missed_cleavages: int = 0) -> int:
    """Digest all proteins and return the total number of tryptic peptides."""
    digestor = oms.ProteaseDigestion()
    digestor.setEnzyme(enzyme)
    digestor.setMissedCleavages(missed_cleavages)
    total = 0
    for entry in entries:
        peptides = []
        digestor.digest(oms.AASequence.fromString(entry.sequence), peptides)
        total += len(peptides)
    return total


def find_duplicates(entries: List[oms.FASTAEntry]) -> List[str]:
    """Find duplicate accession identifiers."""
    seen: Counter = Counter()
    for entry in entries:
        seen[entry.identifier.split()[0]] += 1
    return [acc for acc, count in seen.items() if count > 1]


def compute_statistics(
    input_path: str,
    enzyme: Optional[str] = None,
    missed_cleavages: int = 0,
) -> dict:
    """Compute comprehensive statistics for a FASTA file.

    Returns a dict with protein_count, length_stats, aa_frequency,
    tryptic_peptide_count, and duplicate_accessions.
    """
    entries = load_fasta(input_path)
    stats: dict = {
        "protein_count": len(entries),
        "length_stats": compute_length_stats(entries),
        "aa_frequency": compute_aa_frequency(entries),
        "duplicate_accessions": find_duplicates(entries),
    }
    if enzyme:
        stats["tryptic_peptide_count"] = count_tryptic_peptides(entries, enzyme, missed_cleavages)
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Report statistics for a FASTA database."
    )
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--enzyme", default=None, help="Enzyme for digestion (e.g. Trypsin)")
    parser.add_argument("--missed-cleavages", type=int, default=0, help="Missed cleavages (default: 0)")
    parser.add_argument("--output", default=None, help="Output JSON file (default: stdout)")
    args = parser.parse_args()

    stats = compute_statistics(args.input, args.enzyme, args.missed_cleavages)
    output = json.dumps(stats, indent=2)

    if args.output:
        with open(args.output, "w") as fh:
            fh.write(output + "\n")
        print(f"Statistics written to {args.output}")
    else:
        print(output)


if __name__ == "__main__":
    main()
