"""
Missed Cleavage Analyzer
=========================
Analyze missed cleavage distribution from a peptide list. Useful as a QC metric.

Features
--------
- Count missed cleavages per peptide using enzyme-specific rules
- Generate distribution statistics
- Support for Trypsin, Lys-C, and other common enzymes

Usage
-----
    python missed_cleavage_analyzer.py --input peptides.tsv --enzyme Trypsin --output mc_report.tsv
"""

import argparse
import csv
import json
import sys

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def count_missed_cleavages(sequence: str, enzyme: str = "Trypsin") -> int:
    """Count missed cleavages in a peptide sequence for a given enzyme.

    Parameters
    ----------
    sequence : str
        Peptide sequence.
    enzyme : str
        Enzyme name (e.g., 'Trypsin', 'Lys-C').

    Returns
    -------
    int
        Number of missed cleavages.
    """
    aa_seq = oms.AASequence.fromString(sequence)
    digest = oms.ProteaseDigestion()
    digest.setEnzyme(enzyme)

    # Count internal cleavage sites (K/R for trypsin, not at the end)
    count = digest.missedCleavages(aa_seq)
    return count


def analyze_missed_cleavages(peptides: list, enzyme: str = "Trypsin") -> dict:
    """Analyze missed cleavage distribution for a list of peptides.

    Parameters
    ----------
    peptides : list
        List of peptide sequence strings.
    enzyme : str
        Enzyme name.

    Returns
    -------
    dict
        Dictionary with per-peptide results, distribution, and summary stats.
    """
    results = []
    distribution = {}

    for pep in peptides:
        pep = pep.strip()
        if not pep:
            continue
        mc = count_missed_cleavages(pep, enzyme)
        results.append({"peptide": pep, "missed_cleavages": mc})
        distribution[mc] = distribution.get(mc, 0) + 1

    total = len(results)
    avg_mc = sum(r["missed_cleavages"] for r in results) / total if total > 0 else 0.0
    max_mc = max((r["missed_cleavages"] for r in results), default=0)

    return {
        "enzyme": enzyme,
        "total_peptides": total,
        "average_missed_cleavages": round(avg_mc, 4),
        "max_missed_cleavages": max_mc,
        "distribution": {str(k): v for k, v in sorted(distribution.items())},
        "peptide_results": results,
    }


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(description="Analyze missed cleavage distribution.")
    parser.add_argument("--input", required=True, help="TSV file with 'sequence' column.")
    parser.add_argument("--enzyme", type=str, default="Trypsin", help="Enzyme name (default: Trypsin).")
    parser.add_argument("--output", type=str, help="Output file (.tsv or .json).")
    args = parser.parse_args()

    peptide_list = []
    with open(args.input) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            seq = row.get("sequence", "").strip()
            if seq:
                peptide_list.append(seq)

    analysis = analyze_missed_cleavages(peptide_list, args.enzyme)

    if args.output:
        if args.output.endswith(".json"):
            with open(args.output, "w") as fh:
                json.dump(analysis, fh, indent=2)
        else:
            with open(args.output, "w", newline="") as fh:
                writer = csv.DictWriter(fh, fieldnames=["peptide", "missed_cleavages"], delimiter="\t")
                writer.writeheader()
                writer.writerows(analysis["peptide_results"])
        print(f"Results written to {args.output}")
    else:
        print(f"Enzyme: {analysis['enzyme']}")
        print(f"Total peptides: {analysis['total_peptides']}")
        print(f"Average MC: {analysis['average_missed_cleavages']}")
        print(f"Distribution: {analysis['distribution']}")


if __name__ == "__main__":
    main()
