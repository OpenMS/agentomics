"""
Glycopeptide Mass Calculator
=============================
Calculate glycopeptide masses with glycan compositions.

Built-in glycan residue masses (Da):
- HexNAc = 203.079
- Hex    = 162.053
- Fuc    = 146.058
- NeuAc  = 291.095

Usage
-----
    python glycopeptide_mass_calculator.py --sequence PEPTIDEK --glycan "HexNAc(2)Hex(5)Fuc(1)" --charge 3
    python glycopeptide_mass_calculator.py --sequence PEPTIDEK --glycan "HexNAc(2)Hex(3)" --output masses.tsv
"""

import csv
import re
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit(
        "pyopenms is required. Install it with:  pip install pyopenms"
    )

PROTON = 1.007276

# Built-in glycan residue masses (Da)
GLYCAN_MASSES = {
    "HexNAc": 203.079,
    "Hex": 162.053,
    "Fuc": 146.058,
    "NeuAc": 291.095,
}


def parse_glycan(glycan_str: str) -> dict:
    """Parse a glycan composition string like 'HexNAc(2)Hex(5)Fuc(1)'.

    Parameters
    ----------
    glycan_str:
        Glycan composition string.

    Returns
    -------
    dict
        Mapping of glycan residue name to count.
    """
    pattern = re.compile(r"([A-Za-z]+)\((\d+)\)")
    matches = pattern.findall(glycan_str)
    if not matches:
        raise ValueError(f"Could not parse glycan composition: '{glycan_str}'")
    composition = {}
    for name, count in matches:
        if name not in GLYCAN_MASSES:
            raise ValueError(
                f"Unknown glycan residue '{name}'. Known: {', '.join(GLYCAN_MASSES.keys())}"
            )
        composition[name] = int(count)
    return composition


def glycan_mass(composition: dict) -> float:
    """Calculate total glycan mass from a composition dict.

    Parameters
    ----------
    composition:
        Mapping of glycan residue name to count.

    Returns
    -------
    float
        Total glycan mass in Da.
    """
    total = 0.0
    for name, count in composition.items():
        total += GLYCAN_MASSES[name] * count
    return total


def glycopeptide_mass(
    sequence: str,
    glycan_str: str,
    charge: int = 1,
) -> dict:
    """Calculate glycopeptide mass.

    Parameters
    ----------
    sequence:
        Peptide amino acid sequence.
    glycan_str:
        Glycan composition string, e.g. 'HexNAc(2)Hex(5)Fuc(1)'.
    charge:
        Charge state for m/z.

    Returns
    -------
    dict
        Mass information dictionary.
    """
    aa_seq = oms.AASequence.fromString(sequence)
    peptide_mono = aa_seq.getMonoWeight()

    composition = parse_glycan(glycan_str)
    g_mass = glycan_mass(composition)

    total = peptide_mono + g_mass
    mz = (total + charge * PROTON) / charge

    return {
        "sequence": sequence,
        "glycan": glycan_str,
        "peptide_mass": peptide_mono,
        "glycan_mass": g_mass,
        "glycan_composition": composition,
        "total_mass": total,
        "charge": charge,
        "mz": mz,
    }


def write_tsv(results: list, output_path: str) -> None:
    """Write results to a TSV file.

    Parameters
    ----------
    results:
        List of result dicts.
    output_path:
        Output file path.
    """
    fieldnames = ["sequence", "glycan", "peptide_mass", "glycan_mass", "total_mass", "charge", "mz"]
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in results:
            writer.writerow(row)


@click.command(help="Calculate glycopeptide masses with glycan compositions.")
@click.option("--sequence", required=True, help="Peptide amino acid sequence")
@click.option("--glycan", required=True, help='Glycan composition, e.g. "HexNAc(2)Hex(5)Fuc(1)"')
@click.option("--charge", type=int, default=1, help="Charge state (default: 1)")
@click.option("--output", default=None, help="Output TSV file path")
def main(sequence, glycan, charge, output):
    result = glycopeptide_mass(sequence, glycan, charge=charge)

    print(f"Sequence          : {result['sequence']}")
    print(f"Glycan            : {result['glycan']}")
    print(f"Peptide mass      : {result['peptide_mass']:.6f} Da")
    print(f"Glycan mass       : {result['glycan_mass']:.6f} Da")
    print(f"Total mass        : {result['total_mass']:.6f} Da")
    print(f"Charge            : {result['charge']}+")
    print(f"m/z               : {result['mz']:.6f}")

    if output:
        write_tsv([result], output)
        print(f"\nResults written to {output}")


if __name__ == "__main__":
    main()
