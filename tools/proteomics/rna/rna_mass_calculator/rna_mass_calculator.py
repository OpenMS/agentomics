"""
RNA Mass Calculator
===================
Calculate mass, formula, and isotope patterns for RNA sequences.

Supports standard RNA nucleotides (A, C, G, U). Uses pyopenms NASequence
when available, otherwise falls back to manual calculation using monoisotopic
nucleotide residue masses.

Usage
-----
    python rna_mass_calculator.py --sequence AAUGC --charge 2
    python rna_mass_calculator.py --sequence AAUGCAAUGG --charge 3 --output mass.json
"""

import json
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

PROTON = 1.007276

# Monoisotopic residue masses for RNA nucleotides (internal residue, losing H2O)
# These are the masses of the nucleotide monophosphate residues in an RNA chain.
NUCLEOTIDE_RESIDUE_MASSES = {
    "A": 329.05252,
    "C": 305.04188,
    "G": 345.04744,
    "U": 306.02530,
}

# Water mass added once for the full-length RNA (terminal groups)
WATER_MASS = 18.01056


def _has_na_sequence() -> bool:
    """Check if pyopenms has NASequence support."""
    return hasattr(oms, "NASequence")


def calculate_rna_mass(sequence: str, charge: int = 1) -> dict:
    """Calculate monoisotopic mass and m/z for an RNA sequence.

    Parameters
    ----------
    sequence:
        RNA sequence string using A, C, G, U characters.
    charge:
        Charge state for m/z calculation.

    Returns
    -------
    dict
        Dictionary with sequence, charge, monoisotopic_mass, mz, and formula.
    """
    sequence = sequence.upper().strip()
    for ch in sequence:
        if ch not in NUCLEOTIDE_RESIDUE_MASSES:
            raise ValueError(f"Invalid RNA nucleotide: '{ch}'. Must be one of A, C, G, U.")

    if _has_na_sequence():
        na_seq = oms.NASequence.fromString(sequence)
        mono_mass = na_seq.getMonoWeight()
        formula_str = str(na_seq.getFormula())
    else:
        mono_mass = sum(NUCLEOTIDE_RESIDUE_MASSES[nt] for nt in sequence) + WATER_MASS
        formula_str = _manual_formula(sequence)

    mz = (mono_mass + charge * PROTON) / charge

    return {
        "sequence": sequence,
        "charge": charge,
        "monoisotopic_mass": mono_mass,
        "mz": mz,
        "formula": formula_str,
    }


def _manual_formula(sequence: str) -> str:
    """Compute the molecular formula for an RNA sequence manually.

    Each nucleotide residue contributes its elemental composition. The full
    sequence adds one water molecule for the terminal groups.
    """
    # Elemental compositions of each nucleotide residue (monophosphate, internal)
    compositions = {
        "A": {"C": 10, "H": 12, "N": 5, "O": 6, "P": 1},
        "C": {"C": 9, "H": 12, "N": 3, "O": 7, "P": 1},
        "G": {"C": 10, "H": 12, "N": 5, "O": 7, "P": 1},
        "U": {"C": 9, "H": 11, "N": 2, "O": 8, "P": 1},
    }
    total = {"C": 0, "H": 0, "N": 0, "O": 0, "P": 0}
    for nt in sequence:
        for elem, count in compositions[nt].items():
            total[elem] += count
    # Add water for terminal groups
    total["H"] += 2
    total["O"] += 1
    return "".join(f"{elem}{total[elem]}" for elem in ["C", "H", "N", "O", "P"] if total[elem] > 0)


def calculate_isotope_pattern(sequence: str, n_peaks: int = 5) -> list:
    """Calculate the isotope distribution pattern for an RNA sequence.

    Parameters
    ----------
    sequence:
        RNA sequence string.
    n_peaks:
        Number of isotope peaks to return.

    Returns
    -------
    list
        List of (mass_offset, relative_intensity) tuples.
    """
    sequence = sequence.upper().strip()
    mass_info = calculate_rna_mass(sequence, charge=1)
    formula_str = mass_info["formula"]

    ef = oms.EmpiricalFormula(formula_str)
    isotopes = ef.getIsotopeDistribution(oms.CoarseIsotopePatternGenerator(n_peaks))

    pattern = []
    for iso in isotopes.getContainer():
        pattern.append((iso.getMZ(), iso.getIntensity()))

    return pattern


@click.command(help="Calculate mass/formula/isotopes for RNA sequences.")
@click.option("--sequence", required=True, help="RNA sequence (e.g. AAUGCAAUGG)")
@click.option("--charge", type=int, default=1, help="Charge state for m/z (default: 1)")
@click.option("--isotopes", type=int, default=0, help="Number of isotope peaks to show (default: 0 = off)")
@click.option("--output", default=None, help="Output JSON file (optional)")
def main(sequence, charge, isotopes, output):
    result = calculate_rna_mass(sequence, charge)

    if isotopes > 0:
        pattern = calculate_isotope_pattern(sequence, isotopes)
        result["isotope_pattern"] = [{"mass": m, "intensity": i} for m, i in pattern]

    if output:
        with open(output, "w") as fh:
            json.dump(result, fh, indent=2)
        print(f"Results written to {output}")
    else:
        print(f"Sequence          : {result['sequence']}")
        print(f"Charge            : {result['charge']}+")
        print(f"Monoisotopic mass : {result['monoisotopic_mass']:.6f} Da")
        print(f"m/z               : {result['mz']:.6f}")
        print(f"Formula           : {result['formula']}")
        if "isotope_pattern" in result:
            print("\n--- Isotope Pattern ---")
            for peak in result["isotope_pattern"]:
                print(f"  {peak['mass']:.4f}  {peak['intensity']:.6f}")


if __name__ == "__main__":
    main()
