"""
Lipid Species Resolver
======================
From a sum-composition lipid notation (e.g. PC 36:2), enumerate all possible
acyl chain combinations. For each combination, compute the exact mass using
pyopenms EmpiricalFormula.

Usage
-----
    python lipid_species_resolver.py --input lipids.tsv --lipid-class PC --output resolved.tsv
"""

import csv

import click
import pyopenms as oms

# Head-group formulas for common lipid classes (without acyl chains).
# Acyl chains contribute CnH(2n-2*db-1)O for each chain (ester linkage).
# The head-group formula represents the glycerol backbone + head-group + linkage atoms.
HEADGROUP_FORMULAS = {
    "PC": "C10H18NO8P",   # glycerophosphocholine head group
    "PE": "C5H10NO8P",    # glycerophosphoethanolamine head group (adjusted)
    "PS": "C6H10NO10P",   # glycerophosphoserine head group (adjusted)
    "PG": "C6H11O10P",    # glycerophosphoglycerol head group (adjusted)
    "PI": "C9H15O13P",    # glycerophosphoinositol head group (adjusted)
    "PA": "C3H5O8P",      # phosphatidic acid head group (adjusted)
    "TG": "C3H5O3",       # triacylglycerol backbone (adjusted for 3 chains)
    "DG": "C3H5O3",       # diacylglycerol backbone (adjusted for 2 chains)
}

# Number of acyl chains per lipid class
CHAIN_COUNT = {
    "PC": 2, "PE": 2, "PS": 2, "PG": 2, "PI": 2, "PA": 2,
    "TG": 3, "DG": 2,
}

# Minimum and maximum carbon atoms per single acyl chain
MIN_CHAIN_C = 2
MAX_CHAIN_C = 26
# Maximum double bonds per chain (capped at (carbons - 1) / 2 physiologically)
MAX_CHAIN_DB = 6


def enumerate_chain_combinations(
    total_carbons: int,
    total_double_bonds: int,
    num_chains: int = 2,
    min_c: int = MIN_CHAIN_C,
    max_c: int = MAX_CHAIN_C,
    max_db: int = MAX_CHAIN_DB,
) -> list[list[tuple[int, int]]]:
    """Enumerate acyl chain combinations that sum to given totals.

    Parameters
    ----------
    total_carbons:
        Sum of carbons across all chains.
    total_double_bonds:
        Sum of double bonds across all chains.
    num_chains:
        Number of acyl chains (default 2 for diacyl lipids).
    min_c:
        Minimum carbon atoms per chain.
    max_c:
        Maximum carbon atoms per chain.
    max_db:
        Maximum double bonds per single chain.

    Returns
    -------
    list of list of (carbons, double_bonds)
        Each inner list is one combination of chains, sorted by carbons then db.
        Duplicate combinations (differing only in order) are removed.
    """
    results: list[list[tuple[int, int]]] = []

    def _recurse(remaining_c: int, remaining_db: int, chains_left: int, current: list[tuple[int, int]]) -> None:
        if chains_left == 0:
            if remaining_c == 0 and remaining_db == 0:
                results.append(sorted(current))
            return
        # The minimum/maximum for this chain
        lo_c = max(min_c, remaining_c - max_c * (chains_left - 1))
        hi_c = min(max_c, remaining_c - min_c * (chains_left - 1))
        for c in range(lo_c, hi_c + 1):
            lo_db = max(0, remaining_db - max_db * (chains_left - 1))
            hi_db = min(max_db, remaining_db, (c - 1) // 2 if c > 1 else 0)
            for db in range(lo_db, hi_db + 1):
                _recurse(remaining_c - c, remaining_db - db, chains_left - 1, current + [(c, db)])

    _recurse(total_carbons, total_double_bonds, num_chains, [])

    # Deduplicate (sorted tuples remove order-dependent duplicates)
    seen: set[tuple[tuple[int, int], ...]] = set()
    unique = []
    for combo in results:
        key = tuple(combo)
        if key not in seen:
            seen.add(key)
            unique.append(combo)
    return unique


def acyl_chain_formula(carbons: int, double_bonds: int) -> str:
    """Return the molecular formula contribution of one acyl chain (ester-linked).

    An ester-linked acyl chain contributes C(n)H(2n - 1 - 2*db)O to the lipid.

    Parameters
    ----------
    carbons:
        Number of carbon atoms.
    double_bonds:
        Number of double bonds.

    Returns
    -------
    str
        Formula string, e.g. "C16H31O" for 16:0.
    """
    h_count = 2 * carbons - 1 - 2 * double_bonds
    return f"C{carbons}H{h_count}O"


def lipid_exact_mass(lipid_class: str, chains: list[tuple[int, int]]) -> float:
    """Compute the exact monoisotopic mass of a resolved lipid species.

    Parameters
    ----------
    lipid_class:
        Lipid class abbreviation (e.g. "PC").
    chains:
        List of (carbons, double_bonds) per chain.

    Returns
    -------
    float
        Monoisotopic mass in Da.
    """
    headgroup = HEADGROUP_FORMULAS.get(lipid_class, "")
    if not headgroup:
        return 0.0
    ef = oms.EmpiricalFormula(headgroup)
    for c, db in chains:
        ef = ef + oms.EmpiricalFormula(acyl_chain_formula(c, db))
    return ef.getMonoWeight()


def parse_sum_composition(notation: str) -> tuple[str, int, int]:
    """Parse a sum-composition notation like 'PC 36:2'.

    Parameters
    ----------
    notation:
        String like "PC 36:2" or "TG 54:3".

    Returns
    -------
    tuple of (lipid_class, total_carbons, total_double_bonds)
    """
    parts = notation.strip().split()
    if len(parts) != 2:
        raise ValueError(f"Cannot parse lipid notation: {notation!r}")
    lipid_class = parts[0].upper()
    carbon_db = parts[1].split(":")
    if len(carbon_db) != 2:
        raise ValueError(f"Cannot parse carbon:db from: {parts[1]!r}")
    return lipid_class, int(carbon_db[0]), int(carbon_db[1])


def resolve_lipids(
    lipids: list[dict],
    lipid_class_override: str | None = None,
) -> list[dict]:
    """Resolve sum-composition lipids into acyl chain combinations.

    Parameters
    ----------
    lipids:
        List of dicts with a 'lipid' column containing notation like "PC 36:2".
    lipid_class_override:
        If set, override the lipid class parsed from notation.

    Returns
    -------
    list of dict
        One row per resolved species with chain notation and exact mass.
    """
    results = []
    for row in lipids:
        notation = row.get("lipid", "").strip()
        if not notation:
            continue
        try:
            cls, total_c, total_db = parse_sum_composition(notation)
        except ValueError:
            continue
        if lipid_class_override:
            cls = lipid_class_override.upper()
        num_chains = CHAIN_COUNT.get(cls, 2)
        combos = enumerate_chain_combinations(total_c, total_db, num_chains)
        for combo in combos:
            chain_str = "/".join(f"{c}:{db}" for c, db in combo)
            mass = lipid_exact_mass(cls, combo)
            results.append({
                "input_lipid": notation,
                "lipid_class": cls,
                "resolved_species": f"{cls} {chain_str}",
                "chains": chain_str,
                "exact_mass": round(mass, 6),
            })
    return results


def load_lipids(path: str) -> list[dict]:
    """Load lipid table from TSV. Expects a 'lipid' column."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def write_resolved(resolved: list[dict], path: str) -> None:
    """Write resolved species to TSV."""
    if not resolved:
        with open(path, "w") as fh:
            fh.write("# No resolved species\n")
        return
    fieldnames = list(resolved[0].keys())
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(resolved)


@click.command()
@click.option("--input", "input_file", required=True, help="Lipid table (TSV) with 'lipid' column (e.g. 'PC 36:2')")
@click.option("--lipid-class", default=None, help="Override lipid class (e.g. PC)")
@click.option("--output", required=True, help="Output resolved species (TSV)")
def main(input_file, lipid_class, output) -> None:
    """CLI entry point."""
    lipids = load_lipids(input_file)
    resolved = resolve_lipids(lipids, lipid_class_override=lipid_class)
    write_resolved(resolved, output)
    print(f"Resolved {len(resolved)} species from {len(lipids)} lipid(s), written to {output}")


if __name__ == "__main__":
    main()
