"""
Mass Decomposition Tool
=======================
Find molecular formula compositions for a given mass within tolerance.
Enumerates possible formulas with element constraints using pyopenms EmpiricalFormula.

Features:
- Enumerate possible molecular formulas for a given mass
- Configurable element constraints
- Mass tolerance filtering
- TSV output

Usage
-----
    python mass_decomposition_tool.py --mass 180.0634 --tolerance 0.01
    python mass_decomposition_tool.py --mass 180.0634 --tolerance 0.01 --output decompositions.tsv
"""

import csv
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


# Element exact masses (monoisotopic)
ELEMENT_MASSES = {
    "C": 12.000000,
    "H": 1.0078250,
    "N": 14.003074,
    "O": 15.994915,
    "S": 31.972071,
    "P": 30.973762,
}

# Default element constraints: element -> (min, max)
DEFAULT_CONSTRAINTS = {
    "C": (0, 20),
    "H": (0, 50),
    "N": (0, 5),
    "O": (0, 15),
    "S": (0, 3),
}


def decompose_mass(
    target_mass: float,
    tolerance: float = 0.01,
    constraints: dict[str, tuple[int, int]] | None = None,
) -> list[dict]:
    """Find molecular formula compositions for a given mass within tolerance.

    Parameters
    ----------
    target_mass : float
        Target monoisotopic mass in Da.
    tolerance : float
        Mass tolerance in Da.
    constraints : dict or None
        Element constraints as {element: (min_count, max_count)}.
        Defaults to C:0-20, H:0-50, N:0-5, O:0-15, S:0-3.

    Returns
    -------
    list[dict]
        List of dicts with keys: formula, mass, error_da, C, H, N, O, S.
    """
    if constraints is None:
        constraints = dict(DEFAULT_CONSTRAINTS)

    elements = list(constraints.keys())
    results = []

    # Use iterative approach with pruning
    _enumerate_formulas(target_mass, tolerance, elements, constraints, 0, {}, results)

    # Sort by absolute error
    results.sort(key=lambda x: abs(x["error_da"]))
    return results


def _enumerate_formulas(
    target_mass: float,
    tolerance: float,
    elements: list[str],
    constraints: dict[str, tuple[int, int]],
    elem_idx: int,
    current: dict[str, int],
    results: list[dict],
) -> None:
    """Recursively enumerate formulas with pruning.

    Parameters
    ----------
    target_mass : float
        Target monoisotopic mass.
    tolerance : float
        Mass tolerance in Da.
    elements : list[str]
        List of element symbols.
    constraints : dict
        Element constraints.
    elem_idx : int
        Current element index being enumerated.
    current : dict
        Current element counts.
    results : list[dict]
        Accumulator for valid results.
    """
    if elem_idx == len(elements):
        # Check if current formula matches target mass
        formula_str = _build_formula_string(current)
        if not formula_str:
            return
        try:
            ef = oms.EmpiricalFormula(formula_str)
            mass = ef.getMonoWeight()
            error = abs(mass - target_mass)
            if error <= tolerance:
                result = {
                    "formula": formula_str,
                    "mass": round(mass, 6),
                    "error_da": round(mass - target_mass, 6),
                }
                for elem in elements:
                    result[elem] = current.get(elem, 0)
                results.append(result)
        except Exception:
            pass
        return

    elem = elements[elem_idx]
    min_count, max_count = constraints[elem]
    current_mass = sum(ELEMENT_MASSES.get(e, 0) * current.get(e, 0) for e in elements[:elem_idx])

    for count in range(min_count, max_count + 1):
        test_mass = current_mass + ELEMENT_MASSES.get(elem, 0) * count
        if test_mass > target_mass + tolerance:
            break
        current[elem] = count
        _enumerate_formulas(target_mass, tolerance, elements, constraints, elem_idx + 1, current, results)

    if elem in current:
        del current[elem]


def _build_formula_string(element_counts: dict[str, int]) -> str:
    """Build a molecular formula string from element counts.

    Parameters
    ----------
    element_counts : dict
        Element symbol to count mapping.

    Returns
    -------
    str
        Molecular formula string.
    """
    parts = []
    for elem in ["C", "H", "N", "O", "S", "P"]:
        count = element_counts.get(elem, 0)
        if count > 0:
            parts.append(f"{elem}{count}" if count > 1 else elem)
    return "".join(parts)


def write_tsv(results: list[dict], output_path: str) -> None:
    """Write decomposition results to TSV file.

    Parameters
    ----------
    results : list[dict]
        List of result dictionaries.
    output_path : str
        Path to output TSV file.
    """
    if not results:
        with open(output_path, "w") as f:
            f.write("formula\tmass\terror_da\n")
        return

    fieldnames = list(results[0].keys())
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


@click.command()
@click.option("--mass", type=float, required=True, help="Target mass in Da")
@click.option("--tolerance", type=float, default=0.01, help="Mass tolerance in Da (default: 0.01)")
@click.option("--output", default=None, help="Output TSV file path (default: print to stdout)")
def main(mass, tolerance, output):
    results = decompose_mass(mass, tolerance)

    if output:
        write_tsv(results, output)
        print(f"Wrote {len(results)} decompositions to {output}")
    else:
        if results:
            print("formula\tmass\terror_da")
            for r in results:
                print(f"{r['formula']}\t{r['mass']}\t{r['error_da']}")
        else:
            print("No formulas found within tolerance.")


if __name__ == "__main__":
    main()
