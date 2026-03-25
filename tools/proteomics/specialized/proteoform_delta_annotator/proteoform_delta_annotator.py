"""
Proteoform Delta Annotator
===========================
Annotate mass differences between proteoforms with known PTMs from the
pyopenms ModificationsDB.  Given a TSV of proteoform masses, the tool
computes pairwise or reference-based mass deltas and matches them against
a built-in PTM mass table within a user-specified tolerance.

Usage
-----
    python proteoform_delta_annotator.py --input proteoform_masses.tsv \
        --tolerance 0.5 --output annotated.tsv
"""

import csv
import sys
from typing import Dict, List, Optional, Tuple

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def build_ptm_mass_table() -> List[Dict[str, object]]:
    """Build a table of known PTM names and their monoisotopic mass shifts.

    Returns
    -------
    list of dict
        Each dict has ``name`` (str) and ``mass_shift`` (float).
    """
    mod_db = oms.ModificationsDB()
    mod_names: List[bytes] = []
    mod_db.getAllSearchModifications(mod_names)

    table: List[Dict[str, object]] = []
    seen: set = set()
    for name_bytes in mod_names:
        name = name_bytes.decode("utf-8") if isinstance(name_bytes, bytes) else str(name_bytes)
        if name in seen:
            continue
        seen.add(name)
        try:
            mod = mod_db.getModification(name)
            mass_shift = mod.getDiffMonoMass()
            table.append({"name": name, "mass_shift": mass_shift})
        except Exception:
            continue
    return table


def annotate_delta(
    delta: float, ptm_table: List[Dict[str, object]], tolerance: float
) -> List[Dict[str, object]]:
    """Find PTMs whose mass shift matches *delta* within *tolerance*.

    Parameters
    ----------
    delta:
        Observed mass difference in Da.
    ptm_table:
        Table built by :func:`build_ptm_mass_table`.
    tolerance:
        Absolute mass tolerance in Da.

    Returns
    -------
    list of dict
        Matching PTMs with ``name``, ``mass_shift``, and ``error_da``.
    """
    matches: List[Dict[str, object]] = []
    for entry in ptm_table:
        error = abs(delta - entry["mass_shift"])
        if error <= tolerance:
            matches.append({
                "name": entry["name"],
                "mass_shift": entry["mass_shift"],
                "error_da": error,
            })
    matches.sort(key=lambda m: m["error_da"])
    return matches


def annotate_proteoform_deltas(
    masses: List[Tuple[str, float]],
    tolerance: float,
    reference_mass: Optional[float] = None,
) -> List[Dict[str, object]]:
    """Annotate mass deltas for a list of proteoforms.

    If *reference_mass* is given, each proteoform is compared to it.
    Otherwise the first entry is used as the reference.

    Parameters
    ----------
    masses:
        List of ``(proteoform_id, mass)`` tuples.
    tolerance:
        Absolute tolerance in Da for PTM matching.
    reference_mass:
        Optional reference mass.  Defaults to first entry.

    Returns
    -------
    list of dict
        One entry per proteoform with ``id``, ``mass``, ``delta``, and
        ``annotations`` (semicolon-joined PTM names).
    """
    ptm_table = build_ptm_mass_table()
    if reference_mass is None:
        reference_mass = masses[0][1]

    results: List[Dict[str, object]] = []
    for pf_id, mass in masses:
        delta = mass - reference_mass
        matches = annotate_delta(delta, ptm_table, tolerance)
        annotation_str = "; ".join(
            f"{m['name']} ({m['mass_shift']:.4f} Da, err={m['error_da']:.4f})"
            for m in matches
        )
        results.append({
            "proteoform_id": pf_id,
            "mass": mass,
            "delta": delta,
            "annotations": annotation_str if annotation_str else "no match",
        })
    return results


@click.command(help="Annotate mass differences between proteoforms with known PTMs.")
@click.option("--input", "input", required=True, help="Input TSV with 'proteoform_id' and 'mass' columns")
@click.option("--tolerance", type=float, default=0.5, help="Mass tolerance in Da (default: 0.5)")
@click.option("--output", required=True, help="Output annotated TSV")
def main(input, tolerance, output) -> None:
    masses: List[Tuple[str, float]] = []
    with open(input, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pf_id = row.get("proteoform_id", "").strip()
            mass_str = row.get("mass", "").strip()
            if pf_id and mass_str:
                masses.append((pf_id, float(mass_str)))

    if len(masses) < 1:
        sys.exit("Need at least one proteoform in input.")

    results = annotate_proteoform_deltas(masses, tolerance)

    with open(output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["proteoform_id", "mass", "delta", "annotations"])
        for r in results:
            writer.writerow([r["proteoform_id"], r["mass"], f"{r['delta']:.4f}", r["annotations"]])

    print(f"Annotated {len(results)} proteoforms -> {output}")


if __name__ == "__main__":
    main()
