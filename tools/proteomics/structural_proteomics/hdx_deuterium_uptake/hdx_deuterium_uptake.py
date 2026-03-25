"""
HDX Deuterium Uptake Calculator
================================
Calculate deuterium uptake from HDX-MS data: mass shift, fractional uptake,
and back-exchange correction.

Exchangeable amides = sequence length - number of prolines - 2
(N-terminal amide and first residue do not exchange under typical HDX conditions).

Usage
-----
    python hdx_deuterium_uptake.py --peptides peptides.tsv --undeuterated ref.tsv \
        --timepoints 0,10,60 --output uptake.tsv
"""

import csv
import sys
from typing import Dict, List

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


DEUTERIUM_MASS_SHIFT = 1.00628  # Da difference between D and H


def count_exchangeable_amides(sequence: str) -> int:
    """Count the number of exchangeable backbone amide hydrogens.

    Exchangeable amides = len(sequence) - prolines - 2
    (subtract 2 for N-terminal amino group and first residue which don't exchange).

    Parameters
    ----------
    sequence:
        Amino acid sequence string.

    Returns
    -------
    int
        Number of exchangeable amide hydrogens.
    """
    aa = oms.AASequence.fromString(sequence)
    n = aa.size()
    proline_count = sum(1 for i in range(n) if aa.getResidue(i).getOneLetterCode() == "P")
    exchangeable = n - proline_count - 2
    return max(0, exchangeable)


def get_peptide_mass(sequence: str) -> float:
    """Get the monoisotopic mass of a peptide.

    Parameters
    ----------
    sequence:
        Amino acid sequence.

    Returns
    -------
    float
        Monoisotopic mass in Da.
    """
    aa = oms.AASequence.fromString(sequence)
    return aa.getMonoWeight()


def get_molecular_formula(sequence: str) -> str:
    """Get the molecular formula of a peptide using EmpiricalFormula.

    Parameters
    ----------
    sequence:
        Amino acid sequence.

    Returns
    -------
    str
        Molecular formula string.
    """
    aa = oms.AASequence.fromString(sequence)
    formula = aa.getFormula()
    return formula.toString()


def compute_mass_shift(deuterated_mass: float, undeuterated_mass: float) -> float:
    """Compute the mass shift due to deuterium incorporation.

    Parameters
    ----------
    deuterated_mass:
        Observed centroid mass of deuterated peptide.
    undeuterated_mass:
        Centroid mass of undeuterated reference.

    Returns
    -------
    float
        Mass shift in Da.
    """
    return deuterated_mass - undeuterated_mass


def compute_fractional_uptake(
    mass_shift: float,
    max_exchangeable: int,
    back_exchange_fraction: float = 0.0,
) -> float:
    """Compute fractional deuterium uptake.

    Parameters
    ----------
    mass_shift:
        Observed mass shift in Da.
    max_exchangeable:
        Maximum number of exchangeable amide hydrogens.
    back_exchange_fraction:
        Fraction of deuterium lost to back-exchange (0.0 to 1.0).

    Returns
    -------
    float
        Fractional uptake (0.0 to 1.0).
    """
    if max_exchangeable <= 0:
        return 0.0
    corrected_max = max_exchangeable * (1.0 - back_exchange_fraction)
    if corrected_max <= 0:
        return 0.0
    return mass_shift / (corrected_max * DEUTERIUM_MASS_SHIFT)


def compute_uptake_for_peptide(
    sequence: str,
    undeuterated_mass: float,
    timepoint_masses: Dict[str, float],
    back_exchange_fraction: float = 0.0,
) -> Dict[str, object]:
    """Compute deuterium uptake for a single peptide across timepoints.

    Parameters
    ----------
    sequence:
        Peptide sequence.
    undeuterated_mass:
        Undeuterated reference centroid mass.
    timepoint_masses:
        Dict mapping timepoint labels to observed centroid masses.
    back_exchange_fraction:
        Fraction of back-exchange correction.

    Returns
    -------
    dict
        Dictionary with sequence, formula, exchangeable amides, and per-timepoint uptake.
    """
    exchangeable = count_exchangeable_amides(sequence)
    formula = get_molecular_formula(sequence)
    result: Dict[str, object] = {
        "sequence": sequence,
        "formula": formula,
        "exchangeable_amides": exchangeable,
        "undeuterated_mass": undeuterated_mass,
    }

    for tp, obs_mass in sorted(timepoint_masses.items(), key=lambda x: float(x[0])):
        shift = compute_mass_shift(obs_mass, undeuterated_mass)
        frac = compute_fractional_uptake(shift, exchangeable, back_exchange_fraction)
        result[f"mass_shift_t{tp}"] = round(shift, 6)
        result[f"fractional_uptake_t{tp}"] = round(frac, 6)

    return result


def read_peptides(peptides_path: str) -> List[Dict[str, str]]:
    """Read peptides TSV file with columns: sequence, timepoint, centroid_mass."""
    rows = []
    with open(peptides_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def read_undeuterated(ref_path: str) -> Dict[str, float]:
    """Read undeuterated reference TSV with columns: sequence, centroid_mass."""
    ref = {}
    with open(ref_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            ref[row["sequence"]] = float(row["centroid_mass"])
    return ref


def group_by_peptide(
    rows: List[Dict[str, str]],
) -> Dict[str, Dict[str, float]]:
    """Group peptide rows by sequence, collecting timepoint masses.

    Parameters
    ----------
    rows:
        List of dicts with keys: sequence, timepoint, centroid_mass.

    Returns
    -------
    dict
        Mapping of sequence to {timepoint: centroid_mass}.
    """
    grouped: Dict[str, Dict[str, float]] = {}
    for row in rows:
        seq = row["sequence"]
        tp = row["timepoint"]
        mass = float(row["centroid_mass"])
        if seq not in grouped:
            grouped[seq] = {}
        grouped[seq][tp] = mass
    return grouped


def write_output(output_path: str, results: List[Dict[str, object]]) -> None:
    """Write uptake results to TSV."""
    if not results:
        return
    fieldnames = list(results[0].keys())
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


@click.command(help="Calculate deuterium uptake from HDX-MS data.")
@click.option("--peptides", required=True, help="Peptides TSV (sequence, timepoint, centroid_mass)")
@click.option("--undeuterated", required=True, help="Undeuterated reference TSV (sequence, centroid_mass)")
@click.option("--timepoints", default="0,10,60", help="Comma-separated timepoints to process (default: 0,10,60)")
@click.option("--back-exchange", type=float, default=0.0, help="Back-exchange correction fraction (default: 0.0)")
@click.option("--output", required=True, help="Output uptake TSV file")
def main(peptides, undeuterated, timepoints, back_exchange, output):
    ref = read_undeuterated(undeuterated)
    peptide_rows = read_peptides(peptides)
    grouped = group_by_peptide(peptide_rows)
    timepoints_list = [t.strip() for t in timepoints.split(",")]

    results = []
    for seq, tp_masses in grouped.items():
        if seq not in ref:
            continue
        filtered = {tp: m for tp, m in tp_masses.items() if tp in timepoints_list}
        result = compute_uptake_for_peptide(seq, ref[seq], filtered, back_exchange)
        results.append(result)

    write_output(output, results)

    print(f"Processed {len(results)} peptides")
    print(f"Timepoints: {timepoints}")
    print(f"Back-exchange correction: {back_exchange:.2f}")
    print(f"Output written to {output}")


if __name__ == "__main__":
    main()
