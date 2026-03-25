"""
Isotope Label Detector
=======================
Detect 13C/15N-labeled metabolites by pairing unlabeled and labeled features
based on RT proximity and expected mass shift.

Usage
-----
    python isotope_label_detector.py --unlabeled features_ctrl.tsv \\
        --labeled features_13c.tsv --tracer 13C --ppm 5 --output pairs.tsv
"""

import csv
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


# Mass shift per tracer atom
TRACER_MASS_SHIFTS = {
    "13C": 1.003355,   # 13C - 12C
    "15N": 0.997035,   # 15N - 14N
    "2H": 1.006277,    # 2H - 1H
}


def get_element_count(formula: str, element: str) -> int:
    """Get the count of a specific element in a formula.

    Parameters
    ----------
    formula:
        Molecular formula string.
    element:
        Element symbol (e.g. ``"C"``, ``"N"``).

    Returns
    -------
    int: Count of the element.
    """
    ef = oms.EmpiricalFormula(formula)
    composition = ef.getElementalComposition()
    return composition.get(element.encode(), 0)


def compute_expected_mass_shift(formula: str, tracer: str) -> float:
    """Compute the expected mass shift for full labeling of a formula.

    Parameters
    ----------
    formula:
        Molecular formula string.
    tracer:
        Tracer type (``"13C"``, ``"15N"``, ``"2H"``).

    Returns
    -------
    float: Expected mass shift in Daltons.
    """
    element_map = {"13C": "C", "15N": "N", "2H": "H"}
    element = element_map.get(tracer)
    if element is None:
        raise ValueError(f"Unsupported tracer: {tracer}. Supported: 13C, 15N, 2H")

    n_atoms = get_element_count(formula, element)
    shift_per_atom = TRACER_MASS_SHIFTS[tracer]
    return n_atoms * shift_per_atom


def find_labeled_pairs(
    unlabeled_features: list,
    labeled_features: list,
    tracer: str = "13C",
    ppm: float = 5.0,
    rt_tolerance: float = 10.0,
    max_label_atoms: int = 50,
) -> list:
    """Find pairs of unlabeled and labeled features.

    Parameters
    ----------
    unlabeled_features:
        List of dicts with keys: id, mz, rt (and optionally formula).
    labeled_features:
        List of dicts with keys: id, mz, rt.
    tracer:
        Tracer type.
    ppm:
        Mass tolerance in ppm.
    rt_tolerance:
        RT tolerance in seconds.
    max_label_atoms:
        Maximum number of tracer atoms to consider.

    Returns
    -------
    list of paired feature dicts.
    """
    shift_per_atom = TRACER_MASS_SHIFTS[tracer]
    pairs = []

    for uf in unlabeled_features:
        uf_mz = float(uf["mz"])
        uf_rt = float(uf["rt"])

        for lf in labeled_features:
            lf_mz = float(lf["mz"])
            lf_rt = float(lf["rt"])

            # Check RT proximity
            if abs(uf_rt - lf_rt) > rt_tolerance:
                continue

            # Check if mass difference is a multiple of the tracer shift
            mass_diff = lf_mz - uf_mz
            if mass_diff <= 0:
                continue

            n_labels = round(mass_diff / shift_per_atom)
            if n_labels < 1 or n_labels > max_label_atoms:
                continue

            expected_diff = n_labels * shift_per_atom
            ppm_error = abs(mass_diff - expected_diff) / uf_mz * 1e6

            if ppm_error <= ppm:
                pairs.append({
                    "unlabeled_id": uf.get("id", ""),
                    "unlabeled_mz": round(uf_mz, 6),
                    "labeled_id": lf.get("id", ""),
                    "labeled_mz": round(lf_mz, 6),
                    "mass_diff": round(mass_diff, 6),
                    "n_labels": n_labels,
                    "rt_diff": round(abs(uf_rt - lf_rt), 2),
                    "ppm_error": round(ppm_error, 2),
                })

    return pairs


@click.command()
@click.option("--unlabeled", required=True, help="TSV with unlabeled features (id, mz, rt).")
@click.option("--labeled", required=True, help="TSV with labeled features (id, mz, rt).")
@click.option("--tracer", default="13C", help="Tracer type: 13C, 15N, 2H (default: 13C).")
@click.option("--ppm", type=float, default=5.0, help="Mass tolerance in ppm (default: 5).")
@click.option("--rt-tolerance", type=float, default=10.0,
              help="RT tolerance in seconds (default: 10).")
@click.option("--output", required=True, help="Output TSV with paired features.")
def main(unlabeled, labeled, tracer, ppm, rt_tolerance, output) -> None:
    """CLI entry point."""
    def read_features(path):
        with open(path, newline="") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))

    unlabeled_data = read_features(unlabeled)
    labeled_data = read_features(labeled)

    pairs = find_labeled_pairs(unlabeled_data, labeled_data, tracer=tracer,
                               ppm=ppm, rt_tolerance=rt_tolerance)

    fieldnames = ["unlabeled_id", "unlabeled_mz", "labeled_id", "labeled_mz",
                  "mass_diff", "n_labels", "rt_diff", "ppm_error"]
    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(pairs)

    print(f"Found {len(pairs)} labeled pairs, wrote to {output}")


if __name__ == "__main__":
    main()
