"""
Peptide Property Calculator
===========================
Calculate physicochemical properties of peptide sequences using pyopenms.

Features
--------
- Isoelectric point (pI) via Henderson-Hasselbalch iterative bisection
- GRAVY (grand average of hydropathicity, Kyte-Doolittle scale)
- Net charge at any pH
- Instability index (DIWV weight values)
- Amino acid composition
- Molecular weight and formula

Usage
-----
    python peptide_property_calculator.py --sequence PEPTIDEK --ph 7.0
    python peptide_property_calculator.py --sequence PEPTIDEK --output properties.json
    python peptide_property_calculator.py --input peptides.tsv --output properties.tsv
"""

import csv
import json

import click
import pyopenms as oms

# Kyte-Doolittle hydropathicity scale
KYTE_DOOLITTLE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "E": -3.5, "Q": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# pKa values (Lehninger)
PKA_NTERM = 9.69
PKA_CTERM = 2.34
PKA_SIDE = {
    "C": 8.33, "D": 3.65, "E": 4.25, "H": 6.00, "K": 10.53, "R": 12.48, "Y": 10.07,
}
# +1 for basic, -1 for acidic
CHARGE_SIGN = {
    "C": -1, "D": -1, "E": -1, "H": 1, "K": 1, "R": 1, "Y": -1,
}

# Instability index DIWV weights (Guruprasad et al. 1990) - simplified subset
DIWV = {
    ("D", "G"): 1, ("D", "P"): 1, ("E", "S"): 1,
}


def _charge_at_ph(sequence: str, ph: float) -> float:
    """Calculate net charge at a given pH using Henderson-Hasselbalch.

    Parameters
    ----------
    sequence : str
        One-letter amino acid sequence.
    ph : float
        pH value.

    Returns
    -------
    float
        Estimated net charge.
    """
    charge = 0.0
    # N-terminus contributes +1 when protonated
    charge += 1.0 / (1.0 + 10 ** (ph - PKA_NTERM))
    # C-terminus contributes -1 when deprotonated
    charge -= 1.0 / (1.0 + 10 ** (PKA_CTERM - ph))
    for aa in sequence:
        if aa in PKA_SIDE:
            pka = PKA_SIDE[aa]
            sign = CHARGE_SIGN[aa]
            if sign > 0:
                charge += 1.0 / (1.0 + 10 ** (ph - pka))
            else:
                charge -= 1.0 / (1.0 + 10 ** (pka - ph))
    return charge


def calculate_pi(sequence: str, precision: float = 0.01) -> float:
    """Calculate isoelectric point via bisection on Henderson-Hasselbalch charge.

    Parameters
    ----------
    sequence : str
        One-letter amino acid sequence.
    precision : float
        Desired precision for the pI estimate.

    Returns
    -------
    float
        Estimated isoelectric point.
    """
    low, high = 0.0, 14.0
    while (high - low) > precision:
        mid = (low + high) / 2.0
        c = _charge_at_ph(sequence, mid)
        if c > 0:
            low = mid
        else:
            high = mid
    return round((low + high) / 2.0, 2)


def calculate_gravy(sequence: str) -> float:
    """Calculate GRAVY (Grand Average of Hydropathicity) using Kyte-Doolittle.

    Parameters
    ----------
    sequence : str
        One-letter amino acid sequence.

    Returns
    -------
    float
        GRAVY score.
    """
    values = [KYTE_DOOLITTLE.get(aa, 0.0) for aa in sequence]
    if not values:
        return 0.0
    return round(sum(values) / len(values), 4)


def calculate_instability_index(sequence: str) -> float:
    """Calculate the instability index (simplified Guruprasad method).

    Parameters
    ----------
    sequence : str
        One-letter amino acid sequence.

    Returns
    -------
    float
        Instability index value. Values > 40 suggest the protein is unstable.
    """
    if len(sequence) < 2:
        return 0.0
    total = 0.0
    for i in range(len(sequence) - 1):
        dipeptide = (sequence[i], sequence[i + 1])
        total += DIWV.get(dipeptide, 0.0)
    return round((10.0 / len(sequence)) * total, 4)


def amino_acid_composition(sequence: str) -> dict:
    """Return amino acid counts and frequencies.

    Parameters
    ----------
    sequence : str
        One-letter amino acid sequence.

    Returns
    -------
    dict
        Dictionary with 'counts' and 'frequencies' sub-dicts.
    """
    counts = {}
    for aa in sequence:
        counts[aa] = counts.get(aa, 0) + 1
    length = len(sequence)
    frequencies = {aa: round(count / length, 4) for aa, count in counts.items()} if length else {}
    return {"counts": counts, "frequencies": frequencies}


def calculate_properties(sequence: str, ph: float = 7.0) -> dict:
    """Calculate a full set of physicochemical properties for a peptide.

    Parameters
    ----------
    sequence : str
        Amino acid sequence (plain one-letter code or pyopenms bracket notation).
    ph : float
        pH for net charge calculation.

    Returns
    -------
    dict
        Dictionary containing all computed properties.
    """
    aa_seq = oms.AASequence.fromString(sequence)
    plain = aa_seq.toUnmodifiedString()
    mono = aa_seq.getMonoWeight()
    formula = aa_seq.getFormula()

    pi = calculate_pi(plain)
    gravy = calculate_gravy(plain)
    charge = _charge_at_ph(plain, ph)
    instability = calculate_instability_index(plain)
    composition = amino_acid_composition(plain)

    return {
        "sequence": sequence,
        "unmodified_sequence": plain,
        "length": len(plain),
        "monoisotopic_mass": round(mono, 6),
        "formula": formula.toString(),
        "pI": pi,
        "gravy": gravy,
        "charge_at_ph": round(charge, 4),
        "ph": ph,
        "instability_index": instability,
        "amino_acid_composition": composition,
    }


@click.command(help="Calculate peptide physicochemical properties.")
@click.option("--sequence", type=str, default=None, help="Single peptide sequence.")
@click.option("--ph", type=float, default=7.0, help="pH for charge calculation (default: 7.0).")
@click.option("--input", "input", type=str, default=None, help="TSV file with 'sequence' column.")
@click.option("--output", type=str, default=None, help="Output file (.json or .tsv).")
def main(sequence, ph, input, output):
    """CLI entry point."""
    if not sequence and not input:
        raise click.UsageError("Provide --sequence or --input.")

    results = []
    if sequence:
        results.append(calculate_properties(sequence, ph))
    elif input:
        with open(input) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                seq = row.get("sequence", "").strip()
                if seq:
                    results.append(calculate_properties(seq, ph))

    if output:
        if output.endswith(".json"):
            with open(output, "w") as fh:
                json.dump(results if len(results) > 1 else results[0], fh, indent=2)
        else:
            with open(output, "w", newline="") as fh:
                fieldnames = [
                    "sequence", "unmodified_sequence", "length", "monoisotopic_mass",
                    "formula", "pI", "gravy", "charge_at_ph", "ph", "instability_index",
                ]
                writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
                writer.writeheader()
                for r in results:
                    row = {k: r[k] for k in fieldnames}
                    writer.writerow(row)
        print(f"Results written to {output}")
    else:
        for r in results:
            print(json.dumps(r, indent=2))


if __name__ == "__main__":
    main()
