"""
Drug Metabolite Screener
=========================
Predict drug metabolites from Phase I/II reactions and screen mzML files
for matching ions within a given mass tolerance.

Usage
-----
    python drug_metabolite_screener.py --parent-formula C17H14ClN3O \\
        --reactions phase1,phase2 --input run.mzML --ppm 5 --output metabolites.tsv
"""

import csv
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


# Phase I reactions (functionalization)
PHASE1_REACTIONS = {
    "oxidation": {"add": "O", "remove": ""},
    "demethylation": {"add": "", "remove": "CH2"},
    "hydroxylation": {"add": "O", "remove": ""},
    "dehydrogenation": {"add": "", "remove": "H2"},
    "reduction": {"add": "H2", "remove": ""},
}

# Phase II reactions (conjugation)
PHASE2_REACTIONS = {
    "glucuronidation": {"add": "C6H8O6", "remove": ""},
    "sulfation": {"add": "SO3", "remove": ""},
    "glutathione": {"add": "C10H15N3O6S", "remove": ""},
    "acetylation": {"add": "C2H2O", "remove": ""},
    "methylation": {"add": "CH2", "remove": ""},
}


def get_reaction_table(reaction_sets: list) -> dict:
    """Return the combined reaction table for specified phase sets.

    Parameters
    ----------
    reaction_sets:
        List of phase identifiers, e.g. ``["phase1", "phase2"]``.

    Returns
    -------
    dict mapping reaction names to add/remove formula strings.
    """
    table = {}
    for rs in reaction_sets:
        if rs == "phase1":
            table.update(PHASE1_REACTIONS)
        elif rs == "phase2":
            table.update(PHASE2_REACTIONS)
    return table


def predict_metabolites(parent_formula: str, reaction_sets: list) -> list:
    """Predict metabolite formulas by applying Phase I/II reactions to a parent formula.

    Parameters
    ----------
    parent_formula:
        Molecular formula of the parent drug.
    reaction_sets:
        List of phase identifiers.

    Returns
    -------
    list of dicts with keys: reaction, formula, exact_mass, mass_shift
    """
    parent_ef = oms.EmpiricalFormula(parent_formula)
    parent_mass = parent_ef.getMonoWeight()
    reactions = get_reaction_table(reaction_sets)
    results = []

    for name, mods in reactions.items():
        try:
            new_ef = oms.EmpiricalFormula(parent_formula)
            if mods["add"]:
                add_ef = oms.EmpiricalFormula(mods["add"])
                new_ef = oms.EmpiricalFormula(str(new_ef) + str(add_ef))
            if mods["remove"]:
                remove_ef = oms.EmpiricalFormula(mods["remove"])
                # Build formula by adding negative of removed atoms
                new_ef = oms.EmpiricalFormula(str(new_ef) + "(" + str(remove_ef) + ")-1")

            met_mass = new_ef.getMonoWeight()
            results.append({
                "reaction": name,
                "formula": str(new_ef),
                "exact_mass": round(met_mass, 6),
                "mass_shift": round(met_mass - parent_mass, 6),
            })
        except Exception:
            # Skip reactions that produce invalid formulas
            continue

    return results


def screen_mzml(mzml_path: str, target_masses: list, ppm: float = 5.0) -> list:
    """Screen an mzML file for ions matching target masses.

    Parameters
    ----------
    mzml_path:
        Path to the mzML file.
    target_masses:
        List of dicts with at least 'exact_mass' and 'reaction' keys.
    ppm:
        Mass tolerance in parts per million.

    Returns
    -------
    list of dicts with matched features.
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(mzml_path, exp)

    matches = []
    for spectrum in exp:
        if spectrum.getMSLevel() != 1:
            continue
        rt = spectrum.getRT()
        mzs, intensities = spectrum.get_peaks()

        for target in target_masses:
            target_mz = target["exact_mass"]
            tol = target_mz * ppm / 1e6
            for i, mz in enumerate(mzs):
                if abs(mz - target_mz) <= tol:
                    matches.append({
                        "reaction": target["reaction"],
                        "expected_mass": target["exact_mass"],
                        "observed_mz": round(float(mz), 6),
                        "intensity": round(float(intensities[i]), 1),
                        "rt": round(rt, 2),
                        "ppm_error": round((float(mz) - target_mz) / target_mz * 1e6, 2),
                    })
    return matches


@click.command()
@click.option("--parent-formula", required=True, help="Molecular formula of the parent drug.")
@click.option("--reactions", default="phase1,phase2",
              help="Comma-separated reaction sets: phase1, phase2 (default: phase1,phase2).")
@click.option("--input", "input_file", default=None, help="mzML file to screen (optional).")
@click.option("--ppm", type=float, default=5.0, help="Mass tolerance in ppm (default: 5).")
@click.option("--output", required=True, help="Output TSV file.")
def main(parent_formula, reactions, input_file, ppm, output) -> None:
    """CLI entry point."""
    reaction_sets = [r.strip() for r in reactions.split(",")]
    metabolites = predict_metabolites(parent_formula, reaction_sets)

    if input_file:
        matches = screen_mzml(input_file, metabolites, ppm=ppm)
        fieldnames = ["reaction", "expected_mass", "observed_mz", "intensity", "rt", "ppm_error"]
        output_data = matches
    else:
        fieldnames = ["reaction", "formula", "exact_mass", "mass_shift"]
        output_data = metabolites

    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_data)

    print(f"Wrote {len(output_data)} entries to {output}")


if __name__ == "__main__":
    main()
