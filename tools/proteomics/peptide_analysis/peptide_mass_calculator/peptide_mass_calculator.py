"""
Peptide Mass Calculator
=======================
Calculate monoisotopic and average masses for peptide sequences using pyopenms.

Supports:
- Plain amino acid sequences (e.g. "PEPTIDEK")
- Modified sequences in bracket notation (e.g. "PEPTM[147]IDEK")
- Fragment ion mass series (b-ions and y-ions)
- Multiple charge states

Usage
-----
    python peptide_mass_calculator.py --sequence PEPTIDEK
    python peptide_mass_calculator.py --sequence PEPTM[147]IDEK --charge 2
    python peptide_mass_calculator.py --sequence ACDEFGHIK --fragments
"""


import click
import pyopenms as oms

PROTON = 1.007276


def peptide_masses(sequence: str, charge: int = 1) -> dict:
    """Return monoisotopic and average masses for the given peptide sequence.

    Parameters
    ----------
    sequence:
        Amino acid sequence, optionally with bracket-enclosed modifications,
        e.g. ``"PEPTM[147]IDEK"``.
    charge:
        Desired charge state for m/z calculation (default 1).

    Returns
    -------
    dict
        Dictionary with keys ``monoisotopic_mass``, ``average_mass``,
        ``mz_monoisotopic``, and ``mz_average``.
    """
    aa_seq = oms.AASequence.fromString(sequence)
    mono = aa_seq.getMonoWeight()
    avg = aa_seq.getAverageWeight()
    return {
        "sequence": sequence,
        "charge": charge,
        "monoisotopic_mass": mono,
        "average_mass": avg,
        "mz_monoisotopic": (mono + charge * PROTON) / charge,
        "mz_average": (avg + charge * PROTON) / charge,
    }


def fragment_ions(sequence: str) -> dict:
    """Compute singly charged b-ion and y-ion series for a peptide.

    Parameters
    ----------
    sequence:
        Plain or modified amino acid sequence.

    Returns
    -------
    dict
        Dictionary with keys ``b_ions`` and ``y_ions``, each a list of
        ``(index, mass)`` tuples.
    """
    aa_seq = oms.AASequence.fromString(sequence)
    n = aa_seq.size()

    b_ions = []
    for i in range(1, n):
        prefix = aa_seq.getPrefix(i)
        b_ions.append((i, prefix.getMonoWeight(oms.Residue.ResidueType.BIon, 1)))

    y_ions = []
    for i in range(1, n):
        suffix = aa_seq.getSuffix(i)
        y_ions.append((i, suffix.getMonoWeight(oms.Residue.ResidueType.YIon, 1)))

    return {"b_ions": b_ions, "y_ions": y_ions}


@click.command(help="Calculate peptide/fragment masses using pyopenms.")
@click.option("--sequence", required=True, help="Amino acid sequence (e.g. PEPTIDEK or PEPTM[147]IDEK)")
@click.option("--charge", type=int, default=1, help="Charge state for m/z calculation (default: 1)")
@click.option("--fragments", is_flag=True, help="Also compute b-ion and y-ion series")
def main(sequence, charge, fragments):
    info = peptide_masses(sequence, charge)
    print(f"Sequence          : {info['sequence']}")
    print(f"Charge            : {info['charge']}+")
    print(f"Monoisotopic mass : {info['monoisotopic_mass']:.6f} Da")
    print(f"Average mass      : {info['average_mass']:.6f} Da")
    print(f"m/z (mono)        : {info['mz_monoisotopic']:.6f}")
    print(f"m/z (avg)         : {info['mz_average']:.6f}")

    if fragments:
        ions = fragment_ions(sequence)
        print("\n--- b-ions ---")
        for idx, mass in ions["b_ions"]:
            print(f"  b{idx:>2}  {mass:.6f} Da")
        print("\n--- y-ions ---")
        for idx, mass in ions["y_ions"]:
            print(f"  y{idx:>2}  {mass:.6f} Da")


if __name__ == "__main__":
    main()
