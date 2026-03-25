"""
Top-Down Coverage Calculator
=============================
Compute per-residue bond cleavage coverage from fragment ions in top-down
proteomics.  Given a protein sequence and observed fragment ion masses, the
tool generates theoretical b- and y-ion ladders via pyopenms AASequence,
matches observed masses within a user-specified ppm tolerance, and reports
which backbone bonds were covered.

Usage
-----
    python topdown_coverage_calculator.py --sequence PROTEINSEQ \
        --fragments observed.tsv --tolerance 10 --output coverage.tsv
"""

import argparse
import csv
import sys
from typing import Dict, List, Tuple

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

PROTON = 1.007276


def theoretical_fragments(sequence: str) -> Dict[str, List[Tuple[int, float]]]:
    """Generate theoretical b- and y-ion masses for a protein sequence.

    Parameters
    ----------
    sequence:
        Amino acid sequence string.

    Returns
    -------
    dict
        ``"b"`` and ``"y"`` keys each mapping to a list of
        ``(ion_number, monoisotopic_mass)`` tuples.
    """
    aa = oms.AASequence.fromString(sequence)
    n = aa.size()

    b_ions: List[Tuple[int, float]] = []
    y_ions: List[Tuple[int, float]] = []

    for i in range(1, n):
        # b-ion: prefix of length i
        prefix = aa.getPrefix(i)
        b_mass = prefix.getMonoWeight(oms.Residue.ResidueType.BIon, 1)
        b_ions.append((i, b_mass))

        # y-ion: suffix of length i
        suffix = aa.getSuffix(i)
        y_mass = suffix.getMonoWeight(oms.Residue.ResidueType.YIon, 1)
        y_ions.append((i, y_mass))

    return {"b": b_ions, "y": y_ions}


def match_fragments(
    theoretical: Dict[str, List[Tuple[int, float]]],
    observed_masses: List[float],
    tolerance_ppm: float,
) -> Dict[str, List[Dict[str, object]]]:
    """Match observed masses against theoretical fragment ions.

    Parameters
    ----------
    theoretical:
        Output of :func:`theoretical_fragments`.
    observed_masses:
        List of observed monoisotopic masses (singly charged).
    tolerance_ppm:
        Tolerance in parts-per-million.

    Returns
    -------
    dict
        ``"b"`` and ``"y"`` keys mapping to lists of match dicts with
        ``ion_number``, ``theoretical_mass``, ``observed_mass``, ``error_ppm``.
    """
    matches: Dict[str, List[Dict[str, object]]] = {"b": [], "y": []}
    for ion_type in ("b", "y"):
        for ion_num, theo_mass in theoretical[ion_type]:
            for obs_mass in observed_masses:
                error_ppm = abs(obs_mass - theo_mass) / theo_mass * 1e6
                if error_ppm <= tolerance_ppm:
                    matches[ion_type].append({
                        "ion_number": ion_num,
                        "theoretical_mass": theo_mass,
                        "observed_mass": obs_mass,
                        "error_ppm": error_ppm,
                    })
                    break  # take first match per ion
    return matches


def bond_coverage(
    sequence: str, matches: Dict[str, List[Dict[str, object]]]
) -> List[Dict[str, object]]:
    """Compute per-bond cleavage coverage.

    Bond *i* (between residues *i* and *i+1*, 1-indexed) is covered if
    b_i or y_(n-i) was matched.

    Parameters
    ----------
    sequence:
        Protein sequence.
    matches:
        Output of :func:`match_fragments`.

    Returns
    -------
    list of dict
        One entry per bond with ``bond_index``, ``left_residue``,
        ``right_residue``, ``covered``, ``ion_types``.
    """
    n = len(sequence)
    b_matched = {m["ion_number"] for m in matches["b"]}
    y_matched = {m["ion_number"] for m in matches["y"]}

    coverage: List[Dict[str, object]] = []
    for i in range(1, n):
        ions = []
        if i in b_matched:
            ions.append(f"b{i}")
        y_num = n - i
        if y_num in y_matched:
            ions.append(f"y{y_num}")
        coverage.append({
            "bond_index": i,
            "left_residue": sequence[i - 1],
            "right_residue": sequence[i],
            "covered": len(ions) > 0,
            "ion_types": ",".join(ions) if ions else "-",
        })
    return coverage


def coverage_summary(bond_cov: List[Dict[str, object]]) -> Dict[str, object]:
    """Summarize bond coverage.

    Returns
    -------
    dict
        ``total_bonds``, ``covered_bonds``, ``coverage_fraction``.
    """
    total = len(bond_cov)
    covered = sum(1 for b in bond_cov if b["covered"])
    return {
        "total_bonds": total,
        "covered_bonds": covered,
        "coverage_fraction": covered / total if total > 0 else 0.0,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compute per-residue bond cleavage coverage from fragment ions."
    )
    parser.add_argument("--sequence", required=True, help="Protein amino acid sequence")
    parser.add_argument(
        "--fragments", required=True,
        help="TSV with 'mass' column of observed fragment ion masses",
    )
    parser.add_argument(
        "--tolerance", type=float, default=10.0,
        help="Tolerance in ppm (default: 10)",
    )
    parser.add_argument("--output", required=True, help="Output coverage TSV")
    args = parser.parse_args()

    # Read observed masses
    observed: List[float] = []
    with open(args.fragments, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            mass_str = row.get("mass", "").strip()
            if mass_str:
                observed.append(float(mass_str))

    if not observed:
        sys.exit("No observed masses found in fragments file.")

    theo = theoretical_fragments(args.sequence)
    matches = match_fragments(theo, observed, args.tolerance)
    cov = bond_coverage(args.sequence, matches)
    summary = coverage_summary(cov)

    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["bond_index", "left_residue", "right_residue", "covered", "ion_types"])
        for entry in cov:
            writer.writerow([
                entry["bond_index"], entry["left_residue"],
                entry["right_residue"], entry["covered"], entry["ion_types"],
            ])
        writer.writerow([])
        writer.writerow(["metric", "value"])
        writer.writerow(["total_bonds", summary["total_bonds"]])
        writer.writerow(["covered_bonds", summary["covered_bonds"]])
        writer.writerow(["coverage_fraction", f"{summary['coverage_fraction']:.4f}"])

    print(f"Coverage: {summary['covered_bonds']}/{summary['total_bonds']} "
          f"({summary['coverage_fraction']:.1%}) -> {args.output}")


if __name__ == "__main__":
    main()
