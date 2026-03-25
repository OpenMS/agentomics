"""
FASTA Taxonomy Splitter
=======================
Split a multi-organism FASTA file by taxonomy parsed from headers.

Usage
-----
    python fasta_taxonomy_splitter.py --input combined.fasta --pattern "OS=([^=]+) OX=" --output-dir split/
"""

import os
import re
from collections import defaultdict
from typing import Dict, List, Optional

import click
import pyopenms as oms


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


def extract_taxonomy(entry: oms.FASTAEntry, pattern: str) -> Optional[str]:
    """Extract taxonomy string from a FASTA entry header using the given regex pattern.

    Searches both the identifier and description fields.
    """
    header = f"{entry.identifier} {entry.description}"
    match = re.search(pattern, header)
    if match and match.group(1):
        return match.group(1).strip()
    return None


def sanitize_filename(name: str) -> str:
    """Convert a taxonomy name to a safe filename."""
    safe = re.sub(r"[^\w\s-]", "", name)
    safe = re.sub(r"\s+", "_", safe).strip("_")
    return safe[:100] if safe else "unknown"


def split_by_taxonomy(
    input_path: str,
    output_dir: str,
    pattern: str = r"OS=([^=]+)\s+OX=",
) -> dict:
    """Split a FASTA file by taxonomy extracted from headers.

    Returns statistics about the split.
    """
    entries = load_fasta(input_path)
    groups: Dict[str, List[oms.FASTAEntry]] = defaultdict(list)
    unmatched: List[oms.FASTAEntry] = []

    for entry in entries:
        taxonomy = extract_taxonomy(entry, pattern)
        if taxonomy:
            groups[taxonomy].append(entry)
        else:
            unmatched.append(entry)

    os.makedirs(output_dir, exist_ok=True)

    files_written = {}
    for taxonomy, group_entries in groups.items():
        filename = sanitize_filename(taxonomy) + ".fasta"
        filepath = os.path.join(output_dir, filename)
        save_fasta(group_entries, filepath)
        files_written[taxonomy] = {"file": filepath, "count": len(group_entries)}

    if unmatched:
        unmatched_path = os.path.join(output_dir, "unmatched.fasta")
        save_fasta(unmatched, unmatched_path)
        files_written["_unmatched"] = {"file": unmatched_path, "count": len(unmatched)}

    return {
        "total_entries": len(entries),
        "taxonomy_groups": len(groups),
        "unmatched_count": len(unmatched),
        "files": files_written,
    }


@click.command(help="Split a multi-organism FASTA file by taxonomy from headers.")
@click.option("--input", "input", required=True, help="Input FASTA file")
@click.option(
    "--pattern", default=r"OS=([^=]+)\s+OX=",
    help="Regex pattern with one capture group for taxonomy (default: OS=... OX=)",
)
@click.option("--output-dir", required=True, help="Output directory for split files")
def main(input, pattern, output_dir) -> None:
    stats = split_by_taxonomy(input, output_dir, pattern)
    print(f"Total entries: {stats['total_entries']}")
    print(f"Taxonomy groups: {stats['taxonomy_groups']}")
    print(f"Unmatched: {stats['unmatched_count']}")
    for tax, info in stats["files"].items():
        print(f"  {tax}: {info['count']} entries -> {info['file']}")


if __name__ == "__main__":
    main()
