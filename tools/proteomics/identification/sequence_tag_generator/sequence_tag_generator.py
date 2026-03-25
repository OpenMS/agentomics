"""
Sequence Tag Generator
======================
Generate de novo sequence tags from MS2 spectra by computing mass differences
between sorted peaks and matching to amino acid residue masses.

Usage
-----
    python sequence_tag_generator.py --mz-list "200.1,313.2,426.3,539.4" --intensities "100,200,150,300" \\
        --tolerance 0.02 --min-tag-length 3 --output tags.tsv
"""

import csv
import sys
from typing import List, Optional

import click

try:
    import pyopenms as oms  # noqa: F401
except ImportError:
    sys.exit(
        "pyopenms is required. Install it with:  pip install pyopenms"
    )


def get_residue_masses() -> dict:
    """Build a lookup of amino acid single-letter codes to monoisotopic residue masses.

    Returns
    -------
    dict
        Mapping of single-letter amino acid code to residue mass (Da).
    """
    db = oms.ResidueDB()
    residue_masses = {}
    for code in "ACDEFGHIKLMNPQRSTVWY":
        residue = db.getResidue(code)
        residue_masses[code] = residue.getMonoWeight(oms.Residue.ResidueType.Internal)
    return residue_masses


def match_mass_to_residue(
    mass_diff: float,
    residue_masses: dict,
    tolerance: float = 0.02,
) -> List[str]:
    """Match a mass difference to amino acid residues.

    Parameters
    ----------
    mass_diff:
        Mass difference between two peaks (Da).
    residue_masses:
        Dict mapping residue codes to masses.
    tolerance:
        Mass tolerance in Da.

    Returns
    -------
    list
        List of matching single-letter amino acid codes.
    """
    matches = []
    for code, mass in residue_masses.items():
        if abs(mass_diff - mass) <= tolerance:
            matches.append(code)
    return matches


def generate_tags(
    mz_values: List[float],
    intensities: Optional[List[float]] = None,
    tolerance: float = 0.02,
    min_tag_length: int = 3,
) -> List[dict]:
    """Generate sequence tags from a list of m/z values.

    Sorts peaks by m/z, computes pairwise mass differences between consecutive
    peaks, and builds sequence tags by matching differences to amino acid masses.

    Parameters
    ----------
    mz_values:
        List of m/z values from an MS2 spectrum.
    intensities:
        Optional list of corresponding intensities.
    tolerance:
        Mass tolerance in Da for residue matching.
    min_tag_length:
        Minimum tag length to report (number of residues).

    Returns
    -------
    list
        List of tag dicts with tag sequence, start/end m/z, and length.
    """
    if len(mz_values) < 2:
        return []

    residue_masses = get_residue_masses()

    # Sort peaks by m/z
    if intensities:
        paired = sorted(zip(mz_values, intensities), key=lambda x: x[0])
        sorted_mzs = [p[0] for p in paired]
    else:
        sorted_mzs = sorted(mz_values)

    # Build adjacency: for each pair of peaks, find matching residues
    n = len(sorted_mzs)
    tags = []

    # Use dynamic programming: extend tags greedily
    # For each starting peak, try to build the longest tag
    for start in range(n):
        _extend_tag(start, sorted_mzs, residue_masses, tolerance, min_tag_length, tags, "")

    return tags


def _extend_tag(
    idx: int,
    sorted_mzs: List[float],
    residue_masses: dict,
    tolerance: float,
    min_tag_length: int,
    tags: List[dict],
    current_tag: str,
) -> None:
    """Recursively extend a sequence tag from a given peak index.

    Parameters
    ----------
    idx:
        Current peak index.
    sorted_mzs:
        Sorted list of m/z values.
    residue_masses:
        Residue mass lookup.
    tolerance:
        Mass tolerance.
    min_tag_length:
        Minimum tag length to report.
    tags:
        Accumulator list for found tags.
    current_tag:
        Current tag string being built.
    """
    extended = False
    for next_idx in range(idx + 1, len(sorted_mzs)):
        diff = sorted_mzs[next_idx] - sorted_mzs[idx]
        # Skip if diff is too large to be any amino acid
        if diff > 250:
            break
        matches = match_mass_to_residue(diff, residue_masses, tolerance)
        if matches:
            extended = True
            for aa in matches:
                new_tag = current_tag + aa
                _extend_tag(next_idx, sorted_mzs, residue_masses, tolerance, min_tag_length, tags, new_tag)

    if not extended and len(current_tag) >= min_tag_length:
        tags.append({
            "tag": current_tag,
            "length": len(current_tag),
            "end_mz": round(sorted_mzs[idx], 6),
        })


def write_tsv(tags: List[dict], output_path: str) -> None:
    """Write tags to TSV.

    Parameters
    ----------
    tags:
        List of tag dicts.
    output_path:
        Output file path.
    """
    fieldnames = ["tag", "length", "end_mz"]
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in tags:
            writer.writerow(row)


@click.command(help="Generate de novo sequence tags from MS2 spectra.")
@click.option("--mz-list", required=True, help="Comma-separated list of m/z values")
@click.option("--intensities", default=None, help="Comma-separated list of intensities (optional)")
@click.option("--tolerance", type=float, default=0.02, help="Mass tolerance in Da (default: 0.02)")
@click.option("--min-tag-length", type=int, default=3, help="Minimum tag length (default: 3)")
@click.option("--output", default=None, help="Output TSV file path")
def main(mz_list, intensities, tolerance, min_tag_length, output):
    mz_values = [float(x.strip()) for x in mz_list.split(",")]
    intensities_list = None
    if intensities:
        intensities_list = [float(x.strip()) for x in intensities.split(",")]

    tags = generate_tags(mz_values, intensities_list, tolerance=tolerance, min_tag_length=min_tag_length)

    print(f"Found {len(tags)} sequence tags (min length {min_tag_length})")
    for t in tags[:20]:
        print(f"  {t['tag']} (length {t['length']}, end m/z {t['end_mz']})")
    if len(tags) > 20:
        print(f"  ... and {len(tags) - 20} more")

    if output:
        write_tsv(tags, output)
        print(f"Results written to {output}")


if __name__ == "__main__":
    main()
