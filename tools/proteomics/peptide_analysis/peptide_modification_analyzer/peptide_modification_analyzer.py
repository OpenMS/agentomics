"""
Peptide Modification Analyzer
===============================
Parse modified peptide sequences and output residue-by-residue mass breakdown.

Features
--------
- Parse pyopenms bracket notation for modifications
- Report per-residue monoisotopic mass
- Show modification delta mass per position
- Calculate total and m/z masses

Usage
-----
    python peptide_modification_analyzer.py --sequence "PEPTM(Oxidation)IDE" --charge 2
    python peptide_modification_analyzer.py --sequence "PEPTM(Oxidation)IDE" --output breakdown.tsv
"""

import csv
import json
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

PROTON = 1.007276


def analyze_modification(sequence: str, charge: int = 1) -> dict:
    """Analyze a modified peptide and return residue-by-residue mass breakdown.

    Parameters
    ----------
    sequence : str
        Modified peptide sequence in pyopenms notation (e.g., 'PEPTM(Oxidation)IDE').
    charge : int
        Charge state for m/z calculation.

    Returns
    -------
    dict
        Dictionary with overall mass info and per-residue breakdown.
    """
    aa_seq = oms.AASequence.fromString(sequence)
    total_mono = aa_seq.getMonoWeight()
    mz = (total_mono + charge * PROTON) / charge

    residues = []
    for i in range(aa_seq.size()):
        residue = aa_seq.getResidue(i)
        one_letter = residue.getOneLetterCode()
        mono_mass = residue.getMonoWeight()

        # Check for modification
        mod_name = ""
        mod_mass = 0.0
        if aa_seq.isModified(i):
            mod_name = residue.getModificationName()
            # Unmodified residue mass from ResidueDB
            unmod_residue = oms.ResidueDB().getResidue(one_letter)
            unmod_mass = unmod_residue.getMonoWeight()
            mod_mass = round(mono_mass - unmod_mass, 6)

        residues.append({
            "position": i + 1,
            "residue": one_letter,
            "monoisotopic_mass": round(mono_mass, 6),
            "modification": mod_name,
            "modification_delta_mass": mod_mass,
        })

    return {
        "sequence": sequence,
        "unmodified_sequence": aa_seq.toUnmodifiedString(),
        "length": aa_seq.size(),
        "total_monoisotopic_mass": round(total_mono, 6),
        "charge": charge,
        "mz": round(mz, 6),
        "residue_breakdown": residues,
    }


@click.command(help="Analyze modified peptide residue-by-residue mass breakdown.")
@click.option("--sequence", required=True, help="Modified peptide sequence.")
@click.option("--charge", type=int, default=1, help="Charge state (default: 1).")
@click.option("--output", type=str, default=None, help="Output file (.tsv or .json).")
def main(sequence, charge, output):
    """CLI entry point."""
    result = analyze_modification(sequence, charge)

    if output:
        if output.endswith(".json"):
            with open(output, "w") as fh:
                json.dump(result, fh, indent=2)
        else:
            with open(output, "w", newline="") as fh:
                fieldnames = ["position", "residue", "monoisotopic_mass", "modification",
                              "modification_delta_mass"]
                writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
                writer.writeheader()
                writer.writerows(result["residue_breakdown"])
        print(f"Results written to {output}")
    else:
        print(f"Sequence: {result['sequence']}")
        print(f"Total mass: {result['total_monoisotopic_mass']}")
        print(f"m/z (z={result['charge']}): {result['mz']}")
        print()
        for r in result["residue_breakdown"]:
            mod_str = f" [{r['modification']} +{r['modification_delta_mass']}]" if r["modification"] else ""
            print(f"  {r['position']}\t{r['residue']}\t{r['monoisotopic_mass']}{mod_str}")


if __name__ == "__main__":
    main()
