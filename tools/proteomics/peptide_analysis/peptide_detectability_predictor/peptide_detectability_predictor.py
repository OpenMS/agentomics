"""
Peptide Detectability Predictor
================================
Predict peptide detectability from physicochemical heuristics.

Features
--------
- Digest proteins from FASTA with specified enzyme
- Score peptides based on length, hydrophobicity, charge, and mass
- Rank peptides by predicted detectability

Usage
-----
    python peptide_detectability_predictor.py --input proteins.fasta --enzyme Trypsin --output detectability.tsv
    python peptide_detectability_predictor.py --sequence PEPTIDEK
"""

import csv
import json
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

# Kyte-Doolittle scale for hydrophobicity
KYTE_DOOLITTLE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "E": -3.5, "Q": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}


def calculate_detectability_score(sequence: str) -> dict:
    """Calculate a heuristic detectability score for a peptide.

    Parameters
    ----------
    sequence : str
        Peptide sequence.

    Returns
    -------
    dict
        Dictionary with individual feature scores and overall detectability score.
    """
    aa_seq = oms.AASequence.fromString(sequence)
    plain = aa_seq.toUnmodifiedString()
    length = len(plain)
    mono_mass = aa_seq.getMonoWeight()

    # Feature 1: Length score (optimal 7-25 residues)
    if 7 <= length <= 25:
        length_score = 1.0
    elif length < 7:
        length_score = length / 7.0
    else:
        length_score = max(0.0, 1.0 - (length - 25) / 25.0)

    # Feature 2: Hydrophobicity score (moderate GRAVY preferred)
    gravy = sum(KYTE_DOOLITTLE.get(aa, 0.0) for aa in plain) / length if length > 0 else 0.0
    # Optimal GRAVY around -0.5 to 0.5
    hydro_score = max(0.0, 1.0 - abs(gravy) / 4.0)

    # Feature 3: Mass range score (optimal 800-2500 Da)
    if 800 <= mono_mass <= 2500:
        mass_score = 1.0
    elif mono_mass < 800:
        mass_score = mono_mass / 800.0
    else:
        mass_score = max(0.0, 1.0 - (mono_mass - 2500) / 2500.0)

    # Feature 4: No problematic residues (M, W prone to modification/loss)
    problem_count = sum(1 for aa in plain if aa in "MW")
    problem_score = max(0.0, 1.0 - problem_count * 0.2)

    # Feature 5: Contains at least one basic residue (good for ionization)
    basic_count = sum(1 for aa in plain if aa in "KRH")
    basic_score = min(1.0, basic_count * 0.5)

    overall = round((length_score + hydro_score + mass_score + problem_score + basic_score) / 5.0, 4)

    return {
        "sequence": sequence,
        "length": length,
        "monoisotopic_mass": round(mono_mass, 6),
        "gravy": round(gravy, 4),
        "length_score": round(length_score, 4),
        "hydrophobicity_score": round(hydro_score, 4),
        "mass_score": round(mass_score, 4),
        "problem_residue_score": round(problem_score, 4),
        "basic_residue_score": round(basic_score, 4),
        "detectability_score": overall,
    }


def predict_from_fasta(fasta_path: str, enzyme: str = "Trypsin",
                       missed_cleavages: int = 1, min_length: int = 6,
                       max_length: int = 40) -> list:
    """Digest proteins from FASTA and predict detectability for each peptide.

    Parameters
    ----------
    fasta_path : str
        Path to the FASTA file.
    enzyme : str
        Enzyme for digestion.
    missed_cleavages : int
        Number of allowed missed cleavages.
    min_length : int
        Minimum peptide length.
    max_length : int
        Maximum peptide length.

    Returns
    -------
    list
        List of detectability result dicts, sorted by score descending.
    """
    entries = []
    oms.FASTAFile().load(fasta_path, entries)

    digest = oms.ProteaseDigestion()
    digest.setEnzyme(enzyme)
    digest.setMissedCleavages(missed_cleavages)

    results = []
    seen = set()
    for entry in entries:
        aa_seq = oms.AASequence.fromString(entry.sequence)
        peptides = []
        digest.digest(aa_seq, peptides, min_length, max_length)
        for pep in peptides:
            pep_str = pep.toString()
            if pep_str not in seen:
                seen.add(pep_str)
                score_result = calculate_detectability_score(pep_str)
                score_result["protein"] = entry.identifier
                results.append(score_result)

    results.sort(key=lambda x: x["detectability_score"], reverse=True)
    return results


@click.command(help="Predict peptide detectability from physicochemical heuristics.")
@click.option("--input", "input", type=str, default=None, help="Protein FASTA file.")
@click.option("--sequence", type=str, default=None, help="Single peptide sequence.")
@click.option("--enzyme", type=str, default="Trypsin", help="Enzyme (default: Trypsin).")
@click.option("--missed-cleavages", type=int, default=1, help="Missed cleavages (default: 1).")
@click.option("--output", type=str, default=None, help="Output file (.tsv or .json).")
def main(input, sequence, enzyme, missed_cleavages, output):
    """CLI entry point."""
    if sequence:
        result = calculate_detectability_score(sequence)
        if output:
            with open(output, "w") as fh:
                json.dump(result, fh, indent=2)
        else:
            print(json.dumps(result, indent=2))
    elif input:
        results = predict_from_fasta(input, enzyme, missed_cleavages)
        if output:
            with open(output, "w", newline="") as fh:
                fieldnames = ["sequence", "protein", "length", "monoisotopic_mass", "detectability_score",
                              "length_score", "hydrophobicity_score", "mass_score",
                              "problem_residue_score", "basic_residue_score", "gravy"]
                writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
                writer.writeheader()
                writer.writerows(results)
            print(f"Results written to {output}")
        else:
            for r in results[:20]:
                print(f"{r['sequence']}\t{r['detectability_score']}\t{r['protein']}")
    else:
        raise click.UsageError("Provide --sequence or --input.")


if __name__ == "__main__":
    main()
