"""
RNA Fragment Spectrum Generator
===============================
Generate theoretical RNA fragment spectra including c, y, w, and a-B ion series.

RNA backbone fragmentation follows different rules than peptides. The main
fragment ion types are:
- c ions: 5' fragments (cleavage of the 3'-P-O bond)
- y ions: 3' fragments (cleavage of the 5'-P-O bond)
- w ions: 3' fragments with loss of a base
- a-B ions: 5' fragments with loss of a base

Usage
-----
    python rna_fragment_spectrum_generator.py --sequence AAUGC --charge 2
    python rna_fragment_spectrum_generator.py --sequence AAUGC --charge 1 --output fragments.tsv
"""

import csv

import click
import pyopenms as oms  # noqa: F401

PROTON = 1.007276

# Monoisotopic residue masses for nucleotide residues (internal, losing water between residues)
NUCLEOTIDE_RESIDUE_MASSES = {
    "A": 329.05252,
    "C": 305.04188,
    "G": 345.04744,
    "U": 306.02530,
}

WATER_MASS = 18.01056

# Base masses (free base, for a-B and w ion calculations)
BASE_MASSES = {
    "A": 135.05450,  # adenine
    "C": 111.04326,  # cytosine
    "G": 151.04942,  # guanine
    "U": 112.02728,  # uracil
}


def _prefix_mass(sequence: str, length: int) -> float:
    """Calculate the neutral mass of an RNA prefix of given length."""
    return sum(NUCLEOTIDE_RESIDUE_MASSES[sequence[i]] for i in range(length)) + WATER_MASS


def _suffix_mass(sequence: str, length: int) -> float:
    """Calculate the neutral mass of an RNA suffix of given length."""
    n = len(sequence)
    return sum(NUCLEOTIDE_RESIDUE_MASSES[sequence[n - length + i]] for i in range(length)) + WATER_MASS


def generate_c_ions(sequence: str, charge: int = 1) -> list:
    """Generate c-ion series (5' fragments).

    c ions result from cleavage of the 3'-P-O5' bond, retaining the 5' portion
    with a cyclic phosphate.

    Parameters
    ----------
    sequence:
        RNA sequence.
    charge:
        Charge state.

    Returns
    -------
    list
        List of (ion_label, mz) tuples.
    """
    ions = []
    for i in range(1, len(sequence)):
        # c ion = prefix mass + cyclic phosphate (HPO3 = 79.966)
        mass = _prefix_mass(sequence, i) + 79.96633
        mz = (mass + charge * PROTON) / charge
        ions.append((f"c{i}", mz))
    return ions


def generate_y_ions(sequence: str, charge: int = 1) -> list:
    """Generate y-ion series (3' fragments).

    y ions are the complementary 3' fragments from c-ion cleavage.

    Parameters
    ----------
    sequence:
        RNA sequence.
    charge:
        Charge state.

    Returns
    -------
    list
        List of (ion_label, mz) tuples.
    """
    ions = []
    for i in range(1, len(sequence)):
        mass = _suffix_mass(sequence, i)
        mz = (mass + charge * PROTON) / charge
        ions.append((f"y{i}", mz))
    return ions


def generate_w_ions(sequence: str, charge: int = 1) -> list:
    """Generate w-ion series (3' fragments with base loss).

    w ions = y ions with loss of the 3'-terminal base and water.

    Parameters
    ----------
    sequence:
        RNA sequence.
    charge:
        Charge state.

    Returns
    -------
    list
        List of (ion_label, mz) tuples.
    """
    ions = []
    n = len(sequence)
    for i in range(2, len(sequence)):
        # w ion from position: suffix of length i, lose the 5'-most base of that suffix
        suffix_start = n - i
        base_nt = sequence[suffix_start]
        mass = _suffix_mass(sequence, i) - BASE_MASSES[base_nt] - WATER_MASS
        mz = (mass + charge * PROTON) / charge
        ions.append((f"w{i}", mz))
    return ions


def generate_a_minus_b_ions(sequence: str, charge: int = 1) -> list:
    """Generate a-B ion series (5' fragments with base loss).

    a-B ions = prefix losing the 3'-terminal base and water.

    Parameters
    ----------
    sequence:
        RNA sequence.
    charge:
        Charge state.

    Returns
    -------
    list
        List of (ion_label, mz) tuples.
    """
    ions = []
    for i in range(2, len(sequence)):
        base_nt = sequence[i - 1]
        mass = _prefix_mass(sequence, i) - BASE_MASSES[base_nt] - WATER_MASS
        mz = (mass + charge * PROTON) / charge
        ions.append((f"a-B{i}", mz))
    return ions


def generate_all_fragments(sequence: str, charge: int = 1) -> list:
    """Generate all RNA fragment ion types.

    Parameters
    ----------
    sequence:
        RNA sequence (A, C, G, U).
    charge:
        Charge state.

    Returns
    -------
    list
        List of dicts with keys: ion_type, ion_label, mz, charge.
    """
    sequence = sequence.upper().strip()
    for ch in sequence:
        if ch not in NUCLEOTIDE_RESIDUE_MASSES:
            raise ValueError(f"Invalid RNA nucleotide: '{ch}'.")

    results = []
    for ion_type, gen_func in [("c", generate_c_ions), ("y", generate_y_ions),
                                ("w", generate_w_ions), ("a-B", generate_a_minus_b_ions)]:
        for label, mz in gen_func(sequence, charge):
            results.append({
                "ion_type": ion_type,
                "ion_label": label,
                "mz": mz,
                "charge": charge,
            })
    return results


@click.command(help="Generate theoretical RNA fragment spectra.")
@click.option("--sequence", required=True, help="RNA sequence (e.g. AAUGC)")
@click.option("--charge", type=int, default=1, help="Charge state (default: 1)")
@click.option("--output", default=None, help="Output TSV file (optional)")
def main(sequence, charge, output):
    fragments = generate_all_fragments(sequence, charge)

    if output:
        with open(output, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=["ion_type", "ion_label", "mz", "charge"],
                                    delimiter="\t")
            writer.writeheader()
            writer.writerows(fragments)
        print(f"Wrote {len(fragments)} fragment ions to {output}")
    else:
        print(f"Sequence: {sequence.upper()}")
        print(f"Charge: {charge}+")
        print(f"\n{'Ion':<10} {'Type':<6} {'m/z':>14}")
        print("-" * 32)
        for f in fragments:
            print(f"{f['ion_label']:<10} {f['ion_type']:<6} {f['mz']:>14.4f}")


if __name__ == "__main__":
    main()
