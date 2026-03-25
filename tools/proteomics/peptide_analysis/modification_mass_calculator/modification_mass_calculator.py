"""
Modification Mass Calculator
=============================
Query Unimod modifications by name or mass shift. Compute modified peptide masses.

Features
--------
- Search modifications by name in the Unimod database
- List all available modifications
- Calculate mass of modified peptides
- Report delta mass for modifications

Usage
-----
    python modification_mass_calculator.py --search-mod Phospho
    python modification_mass_calculator.py --list-mods
    python modification_mass_calculator.py --sequence PEPTIDEK --modifications "Oxidation(M):4"
"""

import csv
import json

import click
import pyopenms as oms

PROTON = 1.007276


def search_modification(name: str) -> list:
    """Search for a modification by name in the ModificationsDB.

    Parameters
    ----------
    name : str
        Modification name to search for (e.g., 'Phospho', 'Oxidation').

    Returns
    -------
    list
        List of dicts with modification details.
    """
    mod_db = oms.ModificationsDB()
    mods = set()
    mod_db.searchModifications(mods, name, "", 0)
    results = []
    seen = set()
    for mod in mods:
        full_id = mod.getFullId()
        if full_id in seen:
            continue
        seen.add(full_id)
        results.append({
            "full_id": full_id,
            "name": mod.getId(),
            "delta_mass": round(mod.getDiffMonoMass(), 6),
            "origin": mod.getOrigin(),
        })
    return results


def list_common_modifications() -> list:
    """List commonly used modifications.

    Returns
    -------
    list
        List of dicts with common modification info.
    """
    common = [
        "Oxidation", "Carbamidomethyl", "Phospho", "Acetyl", "Deamidated",
        "Methyl", "Dimethyl", "Trimethyl", "Sulfo", "Nitro",
    ]
    results = []
    for name in common:
        mods = search_modification(name)
        results.extend(mods)
    return results


def modified_peptide_mass(sequence: str, modifications: str = "", charge: int = 1) -> dict:
    """Calculate the mass of a modified peptide.

    Parameters
    ----------
    sequence : str
        Base peptide sequence.
    modifications : str
        Comma-separated modifications in format 'ModName(Residue):position',
        e.g., 'Oxidation(M):4,Phospho(S):7'.
    charge : int
        Charge state for m/z calculation.

    Returns
    -------
    dict
        Dictionary with mass information.
    """
    if modifications:
        # Build modified sequence string
        seq_list = list(sequence)
        mod_entries = []
        for mod_str in modifications.split(","):
            mod_str = mod_str.strip()
            if ":" in mod_str:
                mod_name, pos_str = mod_str.rsplit(":", 1)
                pos = int(pos_str) - 1  # convert to 0-based
                mod_entries.append((pos, mod_name))

        # Sort by position descending to insert from right to left
        mod_entries.sort(key=lambda x: x[0], reverse=True)
        for pos, mod_name in mod_entries:
            # Extract just the mod name without residue for bracket notation
            base_name = mod_name.split("(")[0] if "(" in mod_name else mod_name
            seq_list[pos] = seq_list[pos] + f"({base_name})"

        modified_seq = "".join(seq_list)
    else:
        modified_seq = sequence

    aa_seq = oms.AASequence.fromString(modified_seq)
    mono = aa_seq.getMonoWeight()
    mz = (mono + charge * PROTON) / charge

    return {
        "sequence": sequence,
        "modified_sequence": aa_seq.toString(),
        "charge": charge,
        "monoisotopic_mass": round(mono, 6),
        "mz": round(mz, 6),
    }


@click.command(help="Query Unimod modifications and compute modified peptide masses.")
@click.option("--search-mod", type=str, default=None, help="Search for a modification by name.")
@click.option("--list-mods", is_flag=True, help="List common modifications.")
@click.option("--sequence", type=str, default=None, help="Peptide sequence for mass calculation.")
@click.option("--modifications", type=str, default="", help="Modifications (e.g., 'Oxidation(M):4').")
@click.option("--charge", type=int, default=1, help="Charge state (default: 1).")
@click.option("--output", type=str, default=None, help="Output file.")
def main(search_mod, list_mods, sequence, modifications, charge, output):
    """CLI entry point."""
    if list_mods:
        results = list_common_modifications()
        if output:
            with open(output, "w", newline="") as fh:
                writer = csv.DictWriter(fh, fieldnames=["full_id", "name", "delta_mass", "origin"], delimiter="\t")
                writer.writeheader()
                writer.writerows(results)
        else:
            for r in results:
                print(f"{r['name']}\t{r['delta_mass']}\t{r['origin']}\t{r['full_id']}")
    elif search_mod:
        results = search_modification(search_mod)
        if output:
            with open(output, "w") as fh:
                json.dump(results, fh, indent=2)
        else:
            for r in results:
                print(f"{r['name']}\t{r['delta_mass']}\t{r['origin']}\t{r['full_id']}")
    elif sequence:
        result = modified_peptide_mass(sequence, modifications, charge)
        if output:
            with open(output, "w") as fh:
                json.dump(result, fh, indent=2)
        else:
            print(json.dumps(result, indent=2))
    else:
        click.echo(click.get_current_context().get_help())


if __name__ == "__main__":
    main()
