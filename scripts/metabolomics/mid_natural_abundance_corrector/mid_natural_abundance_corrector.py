"""
MID Natural Abundance Corrector
=================================
Correct mass isotopomer distributions (MIDs) for natural 13C abundance.
Builds a correction matrix from the theoretical isotope distribution and
solves via least-squares to obtain corrected fractional labeling.

Usage
-----
    python mid_natural_abundance_corrector.py --input isotopologues.tsv \\
        --formula C6H12O6 --tracer 13C --output corrected.tsv
"""

import argparse
import csv
import sys

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

try:
    import numpy as np
except ImportError:
    sys.exit("numpy is required. Install it with:  pip install numpy")


# Natural abundance of 13C
NATURAL_13C_ABUNDANCE = 0.01109


def get_num_tracer_atoms(formula: str, tracer: str) -> int:
    """Get the number of tracer element atoms in a formula.

    Parameters
    ----------
    formula:
        Molecular formula string.
    tracer:
        Tracer element, e.g. ``"13C"`` or ``"15N"``.

    Returns
    -------
    int: Number of tracer atoms.
    """
    ef = oms.EmpiricalFormula(formula)
    composition = ef.getElementalComposition()

    element_map = {"13C": b"C", "15N": b"N", "2H": b"H"}
    element_key = element_map.get(tracer)
    if element_key is None:
        raise ValueError(f"Unsupported tracer: {tracer}. Supported: 13C, 15N, 2H")

    return composition.get(element_key, 0)


def build_correction_matrix(n_atoms: int, natural_abundance: float = NATURAL_13C_ABUNDANCE) -> np.ndarray:
    """Build the natural abundance correction matrix.

    The correction matrix C is (n+1) x (n+1) where n is the number of tracer atoms.
    C[i,j] = probability that j labeled atoms produce a signal at isotopomer i
    due to natural abundance of the heavy isotope.

    Parameters
    ----------
    n_atoms:
        Number of tracer element atoms in the molecule.
    natural_abundance:
        Natural abundance of the heavy isotope (default: 0.01109 for 13C).

    Returns
    -------
    numpy.ndarray of shape (n_atoms+1, n_atoms+1).
    """
    n = n_atoms + 1
    C = np.zeros((n, n))

    for j in range(n):
        # j = number of labeled atoms (tracer-derived)
        # remaining unlabeled atoms that can contribute natural abundance signal
        remaining = n_atoms - j
        for i in range(n):
            # i = observed isotopomer (M+i)
            # k = number of naturally labeled atoms from the remaining unlabeled pool
            k = i - j
            if k < 0 or k > remaining:
                C[i, j] = 0.0
            else:
                from math import comb
                C[i, j] = (comb(remaining, k)
                            * (natural_abundance ** k)
                            * ((1 - natural_abundance) ** (remaining - k)))

    return C


def correct_mid(measured_mid: list, formula: str, tracer: str = "13C") -> list:
    """Correct a measured MID for natural isotope abundance.

    Parameters
    ----------
    measured_mid:
        List of measured fractional abundances (M+0, M+1, ..., M+n).
    formula:
        Molecular formula of the metabolite.
    tracer:
        Tracer isotope identifier.

    Returns
    -------
    list of corrected fractional abundances.
    """
    n_atoms = get_num_tracer_atoms(formula, tracer)

    if tracer == "13C":
        abundance = NATURAL_13C_ABUNDANCE
    elif tracer == "15N":
        abundance = 0.00364
    elif tracer == "2H":
        abundance = 0.000115
    else:
        raise ValueError(f"Unsupported tracer: {tracer}")

    # Pad or truncate measured MID to n+1 elements
    mid = np.array(measured_mid[:n_atoms + 1], dtype=float)
    if len(mid) < n_atoms + 1:
        mid = np.pad(mid, (0, n_atoms + 1 - len(mid)))

    C = build_correction_matrix(n_atoms, natural_abundance=abundance)

    # Solve C @ x = mid via least squares, enforce non-negativity
    corrected, _, _, _ = np.linalg.lstsq(C, mid, rcond=None)
    corrected = np.maximum(corrected, 0)

    # Normalize to sum to 1
    total = corrected.sum()
    if total > 0:
        corrected = corrected / total

    return [round(float(v), 6) for v in corrected]


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Correct mass isotopomer distributions for natural 13C abundance."
    )
    parser.add_argument("--input", required=True,
                        help="TSV with columns: sample, M0, M1, M2, ... (fractional abundances).")
    parser.add_argument("--formula", required=True, help="Molecular formula of the metabolite.")
    parser.add_argument("--tracer", default="13C", help="Tracer isotope (default: 13C).")
    parser.add_argument("--output", required=True, help="Output TSV with corrected MIDs.")
    args = parser.parse_args()

    n_atoms = get_num_tracer_atoms(args.formula, args.tracer)

    rows_out = []
    with open(args.input, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        headers = reader.fieldnames or []
        mid_cols = [h for h in headers if h.startswith("M")]
        for row in reader:
            measured = [float(row[c]) for c in mid_cols]
            corrected = correct_mid(measured, args.formula, args.tracer)
            out_row = {"sample": row.get("sample", "")}
            for i, val in enumerate(corrected):
                out_row[f"M{i}_corrected"] = val
            rows_out.append(out_row)

    fieldnames = ["sample"] + [f"M{i}_corrected" for i in range(n_atoms + 1)]
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows_out)

    print(f"Wrote {len(rows_out)} corrected MIDs to {args.output}")


if __name__ == "__main__":
    main()
