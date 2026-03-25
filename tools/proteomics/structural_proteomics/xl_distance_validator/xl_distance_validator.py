"""
Crosslink Distance Validator
=============================
Validate crosslinks against PDB structure distances.

Parses PDB files manually (ATOM lines for CA atoms) to extract residue coordinates,
then computes Euclidean distances between crosslinked residue pairs.

Usage
-----
    python xl_distance_validator.py --crosslinks links.tsv --pdb structure.pdb --max-distance 30 --output distances.tsv
"""

import csv
import math
import sys
from typing import Dict, List, Optional, Tuple

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def parse_pdb_ca_atoms(pdb_path: str) -> Dict[Tuple[str, int], Tuple[float, float, float]]:
    """Parse a PDB file and extract CA (alpha carbon) atom coordinates.

    Parameters
    ----------
    pdb_path:
        Path to the PDB file.

    Returns
    -------
    dict
        Mapping of (chain_id, residue_number) to (x, y, z) coordinates.
    """
    ca_atoms: Dict[Tuple[str, int], Tuple[float, float, float]] = {}

    with open(pdb_path, "r") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue

            chain_id = line[21].strip()
            try:
                res_num = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            except (ValueError, IndexError):
                continue

            ca_atoms[(chain_id, res_num)] = (x, y, z)

    return ca_atoms


def euclidean_distance(
    p1: Tuple[float, float, float],
    p2: Tuple[float, float, float],
) -> float:
    """Compute Euclidean distance between two 3D points.

    Parameters
    ----------
    p1:
        First point (x, y, z).
    p2:
        Second point (x, y, z).

    Returns
    -------
    float
        Euclidean distance in Angstroms.
    """
    return math.sqrt(
        (p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2 + (p1[2] - p2[2]) ** 2
    )


def validate_sequence(sequence: str) -> bool:
    """Validate a peptide sequence using AASequence.

    Parameters
    ----------
    sequence:
        Peptide sequence string.

    Returns
    -------
    bool
        True if the sequence is parseable.
    """
    try:
        oms.AASequence.fromString(sequence)
        return True
    except Exception:
        return False


def validate_crosslink(
    chain1: str,
    res1: int,
    chain2: str,
    res2: int,
    ca_atoms: Dict[Tuple[str, int], Tuple[float, float, float]],
    max_distance: float = 30.0,
) -> Optional[Dict[str, object]]:
    """Validate a single crosslink against PDB structure.

    Parameters
    ----------
    chain1:
        Chain ID of first residue.
    res1:
        Residue number of first crosslinked site.
    chain2:
        Chain ID of second residue.
    res2:
        Residue number of second crosslinked site.
    ca_atoms:
        CA atom coordinate dictionary from parse_pdb_ca_atoms.
    max_distance:
        Maximum allowed CA-CA distance in Angstroms.

    Returns
    -------
    dict or None
        Validation result dict, or None if residues not found in structure.
    """
    key1 = (chain1, res1)
    key2 = (chain2, res2)

    if key1 not in ca_atoms or key2 not in ca_atoms:
        return {
            "chain1": chain1, "residue1": res1,
            "chain2": chain2, "residue2": res2,
            "distance": None,
            "satisfied": "UNKNOWN",
            "note": "Residue(s) not found in PDB",
        }

    dist = euclidean_distance(ca_atoms[key1], ca_atoms[key2])
    satisfied = "YES" if dist <= max_distance else "NO"

    return {
        "chain1": chain1, "residue1": res1,
        "chain2": chain2, "residue2": res2,
        "distance": round(dist, 2),
        "satisfied": satisfied,
        "note": "",
    }


def validate_crosslinks(
    crosslinks: List[Dict[str, str]],
    ca_atoms: Dict[Tuple[str, int], Tuple[float, float, float]],
    max_distance: float = 30.0,
) -> List[Dict[str, object]]:
    """Validate a list of crosslinks against a PDB structure.

    Parameters
    ----------
    crosslinks:
        List of dicts with keys: peptide1, peptide2, chain1, residue1, chain2, residue2.
    ca_atoms:
        CA atom coordinates.
    max_distance:
        Maximum distance threshold.

    Returns
    -------
    list
        List of validation result dicts.
    """
    results = []
    for xl in crosslinks:
        chain1 = xl.get("chain1", "A")
        res1 = int(xl["residue1"])
        chain2 = xl.get("chain2", "A")
        res2 = int(xl["residue2"])

        result = validate_crosslink(chain1, res1, chain2, res2, ca_atoms, max_distance)
        if result is not None:
            # Add peptide info if available
            result["peptide1"] = xl.get("peptide1", "")
            result["peptide2"] = xl.get("peptide2", "")
            result["valid_peptide1"] = str(validate_sequence(xl.get("peptide1", ""))) if xl.get("peptide1") else ""
            result["valid_peptide2"] = str(validate_sequence(xl.get("peptide2", ""))) if xl.get("peptide2") else ""
            results.append(result)

    return results


def read_crosslinks(crosslinks_path: str) -> List[Dict[str, str]]:
    """Read crosslinks TSV file."""
    rows = []
    with open(crosslinks_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def write_output(output_path: str, results: List[Dict[str, object]]) -> None:
    """Write validation results to TSV."""
    if not results:
        return
    fieldnames = list(results[0].keys())
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


@click.command(help="Validate crosslinks against PDB structure distances.")
@click.option("--crosslinks", required=True, help="Crosslinks TSV file")
@click.option("--pdb", required=True, help="PDB structure file")
@click.option(
    "--max-distance", type=float, default=30.0,
    help="Maximum allowed CA-CA distance in Angstroms (default: 30)",
)
@click.option("--output", required=True, help="Output distances TSV file")
def main(crosslinks, pdb, max_distance, output):
    ca_atoms = parse_pdb_ca_atoms(pdb)
    crosslinks_data = read_crosslinks(crosslinks)
    results = validate_crosslinks(crosslinks_data, ca_atoms, max_distance)
    write_output(output, results)

    n_satisfied = sum(1 for r in results if r["satisfied"] == "YES")
    n_violated = sum(1 for r in results if r["satisfied"] == "NO")
    n_unknown = sum(1 for r in results if r["satisfied"] == "UNKNOWN")
    print(f"Total crosslinks: {len(results)}")
    print(f"  Satisfied (dist <= {max_distance} A): {n_satisfied}")
    print(f"  Violated:  {n_violated}")
    print(f"  Unknown:   {n_unknown}")


if __name__ == "__main__":
    main()
