"""
Search Result Merger
====================
Merge multiple identification TSV files with union or intersection consensus.

Each input TSV must have at least a 'peptide' column. Additional columns
(score, protein, etc.) are preserved. The merger identifies PSMs by a
composite key of (peptide, charge, spectrum) or just (peptide) if other
columns are absent.

Usage
-----
    python search_result_merger.py --inputs engine1.tsv engine2.tsv --method union --output merged.tsv
    python search_result_merger.py --inputs engine1.tsv engine2.tsv --method intersection --output merged.tsv
"""

import argparse
import csv
import sys

try:
    import pyopenms as oms  # noqa: F401
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def read_identification_tsv(filepath: str) -> list:
    """Read an identification TSV file.

    Parameters
    ----------
    filepath:
        Path to TSV with at least a 'peptide' column.

    Returns
    -------
    list
        List of dicts (one per row).
    """
    with open(filepath) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def _make_key(row: dict) -> str:
    """Create a composite key for a PSM row."""
    parts = [row.get("peptide", "")]
    if "charge" in row:
        parts.append(str(row["charge"]))
    if "spectrum" in row:
        parts.append(row["spectrum"])
    return "||".join(parts)


def merge_results(input_files: list, method: str = "union") -> tuple:
    """Merge identification results from multiple search engines.

    Parameters
    ----------
    input_files:
        List of TSV file paths.
    method:
        'union' (all PSMs from any engine) or 'intersection' (only PSMs found in all).

    Returns
    -------
    tuple
        (fieldnames, merged_rows) where merged_rows is a list of dicts.
    """
    method = method.lower()
    if method not in ("union", "intersection"):
        raise ValueError(f"Unknown method: '{method}'. Choose 'union' or 'intersection'.")

    all_results = []
    all_fieldnames = []
    for filepath in input_files:
        rows = read_identification_tsv(filepath)
        all_results.append(rows)
        if rows:
            for key in rows[0].keys():
                if key not in all_fieldnames:
                    all_fieldnames.append(key)

    if "source" not in all_fieldnames:
        all_fieldnames.append("source")
    if "n_engines" not in all_fieldnames:
        all_fieldnames.append("n_engines")

    # Build key -> list of (source_index, row)
    key_to_entries = {}
    for file_idx, rows in enumerate(all_results):
        source = input_files[file_idx]
        for row in rows:
            key = _make_key(row)
            if key not in key_to_entries:
                key_to_entries[key] = []
            row_copy = dict(row)
            row_copy["_source_idx"] = file_idx
            row_copy["_source"] = source
            key_to_entries[key].append(row_copy)

    n_files = len(input_files)
    merged = []

    for key, entries in key_to_entries.items():
        source_indices = set(e["_source_idx"] for e in entries)

        if method == "intersection" and len(source_indices) < n_files:
            continue

        # Use the first entry as the base row, add source info
        base = dict(entries[0])
        base.pop("_source_idx", None)
        base.pop("_source", None)
        sources = sorted(set(e["_source"] for e in entries))
        base["source"] = ";".join(sources)
        base["n_engines"] = str(len(source_indices))
        merged.append(base)

    return all_fieldnames, merged


def main():
    parser = argparse.ArgumentParser(description="Merge multiple identification TSV files.")
    parser.add_argument("--inputs", nargs="+", required=True, help="Input TSV files")
    parser.add_argument("--method", default="union", choices=["union", "intersection"],
                        help="Merge method (default: union)")
    parser.add_argument("--output", required=True, help="Output TSV file")
    args = parser.parse_args()

    fieldnames, merged = merge_results(args.inputs, method=args.method)

    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(merged)

    print(f"Method: {args.method}")
    print(f"Input files: {len(args.inputs)}")
    print(f"Merged PSMs: {len(merged)}")
    print(f"Output written to {args.output}")


if __name__ == "__main__":
    main()
