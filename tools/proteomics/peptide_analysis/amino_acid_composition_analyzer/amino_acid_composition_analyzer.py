"""
Amino Acid Composition Analyzer
=================================
Analyze amino acid frequency for proteins in a FASTA file.

Features
--------
- Per-protein amino acid counts and frequencies
- Summary statistics across all proteins
- Molecular weight distribution
- Output in TSV or JSON format

Usage
-----
    python amino_acid_composition_analyzer.py --input proteins.fasta --output composition.tsv
    python amino_acid_composition_analyzer.py --sequence PEPTIDEK --output composition.json
"""

import csv
import json
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

STANDARD_AAS = "ACDEFGHIKLMNPQRSTVWY"


def analyze_composition(sequence: str) -> dict:
    """Analyze amino acid composition of a sequence.

    Parameters
    ----------
    sequence : str
        Amino acid sequence (plain or pyopenms notation).

    Returns
    -------
    dict
        Dictionary with counts, frequencies, and molecular properties.
    """
    aa_seq = oms.AASequence.fromString(sequence)
    plain = aa_seq.toUnmodifiedString()
    length = len(plain)

    counts = {}
    for aa in STANDARD_AAS:
        counts[aa] = 0
    for aa in plain:
        if aa in counts:
            counts[aa] += 1

    frequencies = {}
    for aa, count in counts.items():
        frequencies[aa] = round(count / length, 4) if length > 0 else 0.0

    mono_mass = aa_seq.getMonoWeight()

    # Group properties
    basic = sum(counts.get(aa, 0) for aa in "KRH")
    acidic = sum(counts.get(aa, 0) for aa in "DE")
    hydrophobic = sum(counts.get(aa, 0) for aa in "AILMFWV")
    polar = sum(counts.get(aa, 0) for aa in "STNQ")
    aromatic = sum(counts.get(aa, 0) for aa in "FWY")

    return {
        "sequence": sequence[:50] + "..." if len(sequence) > 50 else sequence,
        "length": length,
        "monoisotopic_mass": round(mono_mass, 6),
        "counts": counts,
        "frequencies": frequencies,
        "basic_residues": basic,
        "acidic_residues": acidic,
        "hydrophobic_residues": hydrophobic,
        "polar_residues": polar,
        "aromatic_residues": aromatic,
    }


def analyze_fasta(fasta_path: str) -> list:
    """Analyze amino acid composition for all proteins in a FASTA file.

    Parameters
    ----------
    fasta_path : str
        Path to FASTA file.

    Returns
    -------
    list
        List of composition result dicts, one per protein.
    """
    entries = []
    oms.FASTAFile().load(fasta_path, entries)

    results = []
    for entry in entries:
        result = analyze_composition(entry.sequence)
        result["accession"] = entry.identifier
        results.append(result)
    return results


@click.command(help="Analyze amino acid composition for FASTA proteins.")
@click.option("--input", "input", type=str, default=None, help="Protein FASTA file.")
@click.option("--sequence", type=str, default=None, help="Single sequence to analyze.")
@click.option("--output", type=str, default=None, help="Output file (.tsv or .json).")
def main(input, sequence, output):
    """CLI entry point."""
    if not input and not sequence:
        raise click.UsageError("Provide --input or --sequence.")

    if sequence:
        results = [analyze_composition(sequence)]
    else:
        results = analyze_fasta(input)

    if output:
        if output.endswith(".json"):
            with open(output, "w") as fh:
                json.dump(results if len(results) > 1 else results[0], fh, indent=2)
        else:
            with open(output, "w", newline="") as fh:
                base_fields = ["accession", "length", "monoisotopic_mass",
                               "basic_residues", "acidic_residues", "hydrophobic_residues",
                               "polar_residues", "aromatic_residues"]
                aa_fields = [f"count_{aa}" for aa in STANDARD_AAS]
                fieldnames = base_fields + aa_fields
                writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
                writer.writeheader()
                for r in results:
                    row = {k: r.get(k, "") for k in base_fields}
                    for aa in STANDARD_AAS:
                        row[f"count_{aa}"] = r["counts"].get(aa, 0)
                    writer.writerow(row)
        print(f"Results written to {output}")
    else:
        for r in results:
            acc = r.get("accession", r.get("sequence", ""))
            print(f"{acc}\tlength={r['length']}\tmass={r['monoisotopic_mass']}")


if __name__ == "__main__":
    main()
