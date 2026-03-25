"""
Modified Peptide Generator
===========================
Generate all modified peptide variants for given variable and fixed modifications.

Features
--------
- Enumerate all possible modification combinations
- Support variable and fixed modifications
- Limit maximum number of simultaneous modifications
- Output modified sequences with masses

Usage
-----
    python modified_peptide_generator.py --sequence PEPTMIDEK --variable-mods Oxidation --max-mods 2
    python modified_peptide_generator.py --sequence PEPTMIDEK --variable-mods Oxidation,Phospho --output variants.tsv
"""

import csv
import json
from itertools import combinations

import click
import pyopenms as oms

PROTON = 1.007276

# Mapping of common modification names to their applicable residues
MOD_RESIDUE_MAP = {
    "Oxidation": ["M", "W"],
    "Phospho": ["S", "T", "Y"],
    "Carbamidomethyl": ["C"],
    "Acetyl": ["K"],
    "Deamidated": ["N", "Q"],
    "Methyl": ["K", "R"],
    "Dimethyl": ["K", "R"],
}


def find_modifiable_sites(sequence: str, mod_name: str) -> list:
    """Find all positions in a sequence where a modification can be applied.

    Parameters
    ----------
    sequence : str
        Plain amino acid sequence.
    mod_name : str
        Modification name (e.g., 'Oxidation').

    Returns
    -------
    list
        List of (position, residue) tuples (1-based positions).
    """
    applicable_residues = MOD_RESIDUE_MAP.get(mod_name, [])
    if not applicable_residues:
        # Try to get from ModificationsDB
        mod_db = oms.ModificationsDB()
        mod_names_list = []
        mod_db.searchModifications(mod_names_list, mod_name, "", oms.ResidueModification.TermSpecificity.ANYWHERE)
        for mn in mod_names_list:
            mod = mod_db.getModification(mn)
            origin = mod.getOrigin()
            if origin and origin not in applicable_residues:
                applicable_residues.append(origin)

    sites = []
    for i, aa in enumerate(sequence):
        if aa in applicable_residues:
            sites.append((i + 1, aa))  # 1-based position
    return sites


def generate_variants(sequence: str, variable_mods: list, fixed_mods: list = None,
                      max_mods: int = 2, charge: int = 1) -> list:
    """Generate all modified peptide variants.

    Parameters
    ----------
    sequence : str
        Plain amino acid sequence.
    variable_mods : list
        List of variable modification names.
    fixed_mods : list
        List of fixed modification names (applied to all applicable sites).
    max_mods : int
        Maximum number of simultaneous variable modifications.
    charge : int
        Charge state for m/z calculation.

    Returns
    -------
    list
        List of dicts with variant information.
    """
    if fixed_mods is None:
        fixed_mods = []

    # Apply fixed modifications first
    base_mods = {}
    for mod_name in fixed_mods:
        sites = find_modifiable_sites(sequence, mod_name)
        for pos, _residue in sites:
            base_mods[pos] = mod_name

    # Collect all variable modification sites
    var_sites = []
    for mod_name in variable_mods:
        sites = find_modifiable_sites(sequence, mod_name)
        for pos, residue in sites:
            if pos not in base_mods:
                var_sites.append((pos, residue, mod_name))

    variants = []
    # Generate combinations of variable mods (0 to max_mods)
    for n in range(0, min(max_mods, len(var_sites)) + 1):
        for combo in combinations(var_sites, n):
            mod_dict = dict(base_mods)
            for pos, _residue, mod_name in combo:
                mod_dict[pos] = mod_name

            # Build modified sequence string
            seq_chars = list(sequence)
            for pos in sorted(mod_dict.keys(), reverse=True):
                mod_name = mod_dict[pos]
                base_name = mod_name.split("(")[0] if "(" in mod_name else mod_name
                seq_chars[pos - 1] = seq_chars[pos - 1] + f"({base_name})"
            modified_seq = "".join(seq_chars)

            aa_seq = oms.AASequence.fromString(modified_seq)
            mono = aa_seq.getMonoWeight()
            mz = (mono + charge * PROTON) / charge

            mod_descriptions = []
            for pos in sorted(mod_dict.keys()):
                mod_descriptions.append(f"{mod_dict[pos]}@{pos}")

            variants.append({
                "sequence": sequence,
                "modified_sequence": aa_seq.toString(),
                "modifications": ";".join(mod_descriptions) if mod_descriptions else "none",
                "num_modifications": len(mod_dict),
                "monoisotopic_mass": round(mono, 6),
                "mz": round(mz, 6),
                "charge": charge,
            })

    return variants


@click.command(help="Generate modified peptide variants.")
@click.option("--sequence", required=True, help="Peptide sequence.")
@click.option("--variable-mods", type=str, default="", help="Comma-separated variable modification names.")
@click.option("--fixed-mods", type=str, default="", help="Comma-separated fixed modification names.")
@click.option("--max-mods", type=int, default=2, help="Maximum simultaneous variable mods (default: 2).")
@click.option("--charge", type=int, default=1, help="Charge state (default: 1).")
@click.option("--output", type=str, default=None, help="Output file (.tsv or .json).")
def main(sequence, variable_mods, fixed_mods, max_mods, charge, output):
    """CLI entry point."""
    var_mods = [m.strip() for m in variable_mods.split(",") if m.strip()]
    fix_mods = [m.strip() for m in fixed_mods.split(",") if m.strip()]

    variants = generate_variants(sequence, var_mods, fix_mods, max_mods, charge)

    if output:
        if output.endswith(".json"):
            with open(output, "w") as fh:
                json.dump(variants, fh, indent=2)
        else:
            with open(output, "w", newline="") as fh:
                fieldnames = ["sequence", "modified_sequence", "modifications", "num_modifications",
                              "monoisotopic_mass", "mz", "charge"]
                writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
                writer.writeheader()
                writer.writerows(variants)
        print(f"Generated {len(variants)} variants -> {output}")
    else:
        for v in variants:
            print(f"{v['modified_sequence']}\t{v['modifications']}\t{v['monoisotopic_mass']}")


if __name__ == "__main__":
    main()
