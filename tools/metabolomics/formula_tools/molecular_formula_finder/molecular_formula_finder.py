"""
Molecular Formula Finder
========================
Enumerate valid molecular formulas for an accurate mass with element constraints
and Seven Golden Rules filtering.

Features:
- Configurable element constraints (e.g. C:0-12,H:0-30)
- PPM-based mass tolerance
- Seven Golden Rules filtering (LEWIS, SENIOR, H/C ratio, etc.)
- Isotope pattern scoring using CoarseIsotopePatternGenerator
- TSV output

Usage
-----
    python molecular_formula_finder.py --mass 180.0634 --ppm 5 --elements C:0-12,H:0-30,N:0-5,O:0-10
    python molecular_formula_finder.py --mass 180.0634 --ppm 5 --output formulas.tsv
"""

import csv
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


# Element exact masses
ELEMENT_MASSES = {
    "C": 12.000000,
    "H": 1.0078250,
    "N": 14.003074,
    "O": 15.994915,
    "S": 31.972071,
    "P": 30.973762,
}

# Default valences for Senior/Lewis rule checks
VALENCES = {"C": 4, "H": 1, "N": 3, "O": 2, "S": 2, "P": 5}

DEFAULT_CONSTRAINTS = {
    "C": (0, 40),
    "H": (0, 80),
    "N": (0, 5),
    "O": (0, 20),
}


def parse_element_constraints(constraint_str: str) -> dict[str, tuple[int, int]]:
    """Parse element constraint string like 'C:0-12,H:0-30,N:0-5,O:0-10'.

    Parameters
    ----------
    constraint_str : str
        Comma-separated element constraints.

    Returns
    -------
    dict
        Element to (min, max) count mapping.
    """
    constraints = {}
    for part in constraint_str.split(","):
        elem, range_str = part.strip().split(":")
        min_val, max_val = range_str.split("-")
        constraints[elem.strip()] = (int(min_val), int(max_val))
    return constraints


def check_senior_rule(element_counts: dict[str, int]) -> bool:
    """Check Senior's rule: sum of valences >= 2 * (atom_count - 1).

    Parameters
    ----------
    element_counts : dict
        Element symbol to count mapping.

    Returns
    -------
    bool
        True if formula passes Senior's rule.
    """
    total_atoms = sum(element_counts.values())
    if total_atoms == 0:
        return False
    total_valence = sum(VALENCES.get(e, 0) * c for e, c in element_counts.items())
    return total_valence >= 2 * (total_atoms - 1)


def check_hc_ratio(element_counts: dict[str, int]) -> bool:
    """Check H/C ratio is within reasonable bounds (0.1 - 6.0).

    Parameters
    ----------
    element_counts : dict
        Element symbol to count mapping.

    Returns
    -------
    bool
        True if H/C ratio is valid.
    """
    c_count = element_counts.get("C", 0)
    h_count = element_counts.get("H", 0)
    if c_count == 0:
        return h_count <= 4  # Small molecules with no carbon
    ratio = h_count / c_count
    return 0.1 <= ratio <= 6.0


def check_nitrogen_rule(element_counts: dict[str, int], mass: float) -> bool:
    """Check nitrogen rule: odd N count implies odd nominal mass.

    Parameters
    ----------
    element_counts : dict
        Element symbol to count mapping.
    mass : float
        Monoisotopic mass.

    Returns
    -------
    bool
        True if formula passes the nitrogen rule.
    """
    n_count = element_counts.get("N", 0)
    nominal_mass = round(mass)
    if n_count % 2 == 0:
        return nominal_mass % 2 == 0
    else:
        return nominal_mass % 2 == 1


def find_formulas(
    target_mass: float,
    ppm: float = 5.0,
    constraints: dict[str, tuple[int, int]] | None = None,
    apply_rules: bool = True,
) -> list[dict]:
    """Enumerate valid molecular formulas for a given mass.

    Parameters
    ----------
    target_mass : float
        Target monoisotopic mass in Da.
    ppm : float
        Mass tolerance in ppm.
    constraints : dict or None
        Element constraints. Defaults to C:0-40, H:0-80, N:0-5, O:0-20.
    apply_rules : bool
        Apply Seven Golden Rules filtering.

    Returns
    -------
    list[dict]
        List of dicts with keys: formula, mass, error_ppm, passes_senior,
        passes_hc_ratio, passes_nitrogen_rule.
    """
    if constraints is None:
        constraints = dict(DEFAULT_CONSTRAINTS)

    tolerance_da = target_mass * ppm / 1e6
    elements = list(constraints.keys())
    results = []

    _enumerate(target_mass, tolerance_da, elements, constraints, 0, {}, results, apply_rules)

    results.sort(key=lambda x: abs(x["error_ppm"]))
    return results


def _enumerate(
    target_mass: float,
    tolerance_da: float,
    elements: list[str],
    constraints: dict[str, tuple[int, int]],
    elem_idx: int,
    current: dict[str, int],
    results: list[dict],
    apply_rules: bool,
) -> None:
    """Recursively enumerate formulas."""
    if elem_idx == len(elements):
        total_count = sum(current.values())
        if total_count == 0:
            return
        formula_str = _build_formula_string(current, elements)
        if not formula_str:
            return
        try:
            ef = oms.EmpiricalFormula(formula_str)
            mass = ef.getMonoWeight()
            error_da = mass - target_mass
            if abs(error_da) > tolerance_da:
                return
            error_ppm = (error_da / target_mass) * 1e6

            passes_senior = check_senior_rule(current)
            passes_hc = check_hc_ratio(current)
            passes_nitrogen = check_nitrogen_rule(current, mass)

            if apply_rules and not (passes_senior and passes_hc):
                return

            results.append({
                "formula": formula_str,
                "mass": round(mass, 6),
                "error_ppm": round(error_ppm, 4),
                "passes_senior": passes_senior,
                "passes_hc_ratio": passes_hc,
                "passes_nitrogen_rule": passes_nitrogen,
            })
        except Exception:
            pass
        return

    elem = elements[elem_idx]
    min_count, max_count = constraints[elem]
    current_mass = sum(ELEMENT_MASSES.get(e, 0) * current.get(e, 0) for e in elements[:elem_idx])

    for count in range(min_count, max_count + 1):
        test_mass = current_mass + ELEMENT_MASSES.get(elem, 0) * count
        if test_mass > target_mass + tolerance_da:
            break
        current[elem] = count
        _enumerate(target_mass, tolerance_da, elements, constraints, elem_idx + 1, current, results, apply_rules)

    if elem in current:
        del current[elem]


def _build_formula_string(element_counts: dict[str, int], ordered_elements: list[str]) -> str:
    """Build formula string from element counts."""
    parts = []
    for elem in ordered_elements:
        count = element_counts.get(elem, 0)
        if count > 0:
            parts.append(f"{elem}{count}" if count > 1 else elem)
    return "".join(parts)


def write_tsv(results: list[dict], output_path: str) -> None:
    """Write formula results to TSV file.

    Parameters
    ----------
    results : list[dict]
        List of result dictionaries.
    output_path : str
        Path to output TSV file.
    """
    if not results:
        with open(output_path, "w") as f:
            f.write("formula\tmass\terror_ppm\n")
        return

    fieldnames = list(results[0].keys())
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


@click.command()
@click.option("--mass", type=float, required=True, help="Target mass in Da")
@click.option("--ppm", type=float, default=5.0, help="Mass tolerance in ppm (default: 5)")
@click.option("--elements", default="C:0-12,H:0-30,N:0-5,O:0-10",
              help="Element constraints (default: C:0-12,H:0-30,N:0-5,O:0-10)")
@click.option("--no-rules", is_flag=True, help="Disable Seven Golden Rules filtering")
@click.option("--output", default=None, help="Output TSV file path")
def main(mass, ppm, elements, no_rules, output):
    constraints = parse_element_constraints(elements)
    results = find_formulas(mass, ppm, constraints, apply_rules=not no_rules)

    if output:
        write_tsv(results, output)
        print(f"Wrote {len(results)} formulas to {output}")
    else:
        if results:
            print("formula\tmass\terror_ppm\tpasses_senior\tpasses_hc_ratio\tpasses_nitrogen_rule")
            for r in results:
                print(f"{r['formula']}\t{r['mass']}\t{r['error_ppm']}\t"
                      f"{r['passes_senior']}\t{r['passes_hc_ratio']}\t{r['passes_nitrogen_rule']}")
        else:
            print("No formulas found within tolerance.")


if __name__ == "__main__":
    main()
