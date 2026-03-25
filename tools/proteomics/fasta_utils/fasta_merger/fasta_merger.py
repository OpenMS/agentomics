"""
FASTA Merger
============
Merge multiple FASTA files with optional deduplication.

Usage
-----
    python fasta_merger.py --inputs db1.fasta db2.fasta --remove-duplicates --output merged.fasta
"""

import sys
from typing import List

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def load_fasta(input_path: str) -> List[oms.FASTAEntry]:
    """Load entries from a FASTA file."""
    entries = []
    fasta_file = oms.FASTAFile()
    fasta_file.load(input_path, entries)
    return entries


def save_fasta(entries: List[oms.FASTAEntry], output_path: str) -> None:
    """Save entries to a FASTA file."""
    fasta_file = oms.FASTAFile()
    fasta_file.store(output_path, entries)


def deduplicate_by_identifier(entries: List[oms.FASTAEntry]) -> List[oms.FASTAEntry]:
    """Remove entries with duplicate identifiers, keeping the first occurrence."""
    seen = set()
    result = []
    for entry in entries:
        key = entry.identifier.split()[0]
        if key not in seen:
            seen.add(key)
            result.append(entry)
    return result


def deduplicate_by_sequence(entries: List[oms.FASTAEntry]) -> List[oms.FASTAEntry]:
    """Remove entries with duplicate sequences, keeping the first occurrence."""
    seen = set()
    result = []
    for entry in entries:
        if entry.sequence not in seen:
            seen.add(entry.sequence)
            result.append(entry)
    return result


def merge_fasta_files(
    input_paths: List[str],
    output_path: str,
    remove_duplicates: bool = False,
    dedup_by: str = "identifier",
) -> dict:
    """Merge multiple FASTA files into one.

    Parameters
    ----------
    input_paths : list of str
        Paths to input FASTA files.
    output_path : str
        Path to output FASTA file.
    remove_duplicates : bool
        Whether to remove duplicates.
    dedup_by : str
        Deduplication strategy: 'identifier' or 'sequence'.

    Returns
    -------
    dict
        Statistics about the merge.
    """
    all_entries = []
    file_counts = {}
    for path in input_paths:
        entries = load_fasta(path)
        file_counts[path] = len(entries)
        all_entries.extend(entries)

    total_before = len(all_entries)

    if remove_duplicates:
        if dedup_by == "sequence":
            all_entries = deduplicate_by_sequence(all_entries)
        else:
            all_entries = deduplicate_by_identifier(all_entries)

    save_fasta(all_entries, output_path)
    return {
        "file_counts": file_counts,
        "total_before_dedup": total_before,
        "total_output": len(all_entries),
    }


@click.command(help="Merge multiple FASTA files.")
@click.option("--inputs", multiple=True, required=True, help="Input FASTA files")
@click.option("--output", required=True, help="Output merged FASTA file")
@click.option("--remove-duplicates", is_flag=True, help="Remove duplicate entries")
@click.option(
    "--dedup-by", type=click.Choice(["identifier", "sequence"]),
    default="identifier", help="Deduplication criterion (default: identifier)",
)
def main(inputs, output, remove_duplicates, dedup_by) -> None:
    stats = merge_fasta_files(list(inputs), output, remove_duplicates, dedup_by)
    print(f"Merged {stats['total_before_dedup']} entries -> {stats['total_output']} written to {output}")


if __name__ == "__main__":
    main()
