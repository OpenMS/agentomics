"""
FASTA Cleaner
=============
Clean a FASTA database by removing duplicates, fixing headers, filtering by length,
and removing stop codons.

Usage
-----
    python fasta_cleaner.py --input messy.fasta --remove-duplicates --min-length 6 --output clean.fasta
"""

import argparse
import re
import sys
from typing import List, Optional

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


def remove_duplicates(entries: List[oms.FASTAEntry]) -> List[oms.FASTAEntry]:
    """Remove entries with duplicate sequences."""
    seen = set()
    result = []
    for entry in entries:
        if entry.sequence not in seen:
            seen.add(entry.sequence)
            result.append(entry)
    return result


def remove_stop_codons(entries: List[oms.FASTAEntry]) -> List[oms.FASTAEntry]:
    """Remove trailing stop codons (*) from sequences."""
    for entry in entries:
        entry.sequence = entry.sequence.rstrip("*")
    return entries


def fix_headers(entries: List[oms.FASTAEntry]) -> List[oms.FASTAEntry]:
    """Fix common header issues: remove extra whitespace, sanitize control characters."""
    for entry in entries:
        entry.identifier = re.sub(r"\s+", " ", entry.identifier).strip()
        entry.description = re.sub(r"\s+", " ", entry.description).strip()
    return entries


def filter_by_length(
    entries: List[oms.FASTAEntry],
    min_length: Optional[int] = None,
    max_length: Optional[int] = None,
) -> List[oms.FASTAEntry]:
    """Filter entries by sequence length."""
    result = []
    for entry in entries:
        seq_len = len(entry.sequence)
        if min_length is not None and seq_len < min_length:
            continue
        if max_length is not None and seq_len > max_length:
            continue
        result.append(entry)
    return result


def remove_invalid_chars(entries: List[oms.FASTAEntry]) -> List[oms.FASTAEntry]:
    """Remove non-amino-acid characters from sequences."""
    valid = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
    for entry in entries:
        entry.sequence = "".join(c for c in entry.sequence.upper() if c in valid)
    return entries


def clean_fasta(
    input_path: str,
    output_path: str,
    dedup: bool = False,
    min_length: Optional[int] = None,
    max_length: Optional[int] = None,
    strip_stop_codons: bool = False,
    do_fix_headers: bool = False,
    do_remove_invalid: bool = False,
) -> dict:
    """Clean a FASTA file with the specified operations.

    Returns statistics about the cleaning process.
    """
    entries = load_fasta(input_path)
    total_input = len(entries)

    if strip_stop_codons:
        entries = remove_stop_codons(entries)

    if do_remove_invalid:
        entries = remove_invalid_chars(entries)

    if do_fix_headers:
        entries = fix_headers(entries)

    if dedup:
        entries = remove_duplicates(entries)

    if min_length is not None or max_length is not None:
        entries = filter_by_length(entries, min_length, max_length)

    save_fasta(entries, output_path)
    return {"total_input": total_input, "total_output": len(entries)}


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Clean a FASTA database: remove duplicates, fix headers, filter by length."
    )
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--output", required=True, help="Output cleaned FASTA file")
    parser.add_argument("--remove-duplicates", action="store_true", help="Remove duplicate sequences")
    parser.add_argument("--min-length", type=int, default=None, help="Minimum sequence length")
    parser.add_argument("--max-length", type=int, default=None, help="Maximum sequence length")
    parser.add_argument("--remove-stop-codons", action="store_true", help="Remove trailing stop codons (*)")
    parser.add_argument("--fix-headers", action="store_true", help="Fix header whitespace issues")
    parser.add_argument("--remove-invalid-chars", action="store_true", help="Remove non-amino-acid characters")
    args = parser.parse_args()

    stats = clean_fasta(
        args.input,
        args.output,
        dedup=args.remove_duplicates,
        min_length=args.min_length,
        max_length=args.max_length,
        strip_stop_codons=args.remove_stop_codons,
        do_fix_headers=args.fix_headers,
        do_remove_invalid=args.remove_invalid_chars,
    )
    print(f"Cleaned: {stats['total_input']} -> {stats['total_output']} proteins written to {args.output}")


if __name__ == "__main__":
    main()
