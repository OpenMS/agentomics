"""
FASTA Subset Extractor
======================
Extract proteins from a FASTA database by accession list, keyword, or length range.

Usage
-----
    python fasta_subset_extractor.py --input db.fasta --accessions list.txt --output subset.fasta
    python fasta_subset_extractor.py --input db.fasta --keyword "Homo sapiens" --output subset.fasta
    python fasta_subset_extractor.py --input db.fasta --min-length 50 --max-length 500 --output subset.fasta
"""

import sys
from typing import List, Optional

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


def filter_by_accessions(entries: List[oms.FASTAEntry], accessions: set) -> List[oms.FASTAEntry]:
    """Filter FASTA entries by a set of accession identifiers."""
    result = []
    for entry in entries:
        identifier = entry.identifier.split()[0]
        # Also try extracting UniProt-style accession: sp|P12345|NAME
        parts = identifier.split("|")
        accession_variants = {identifier}
        if len(parts) >= 2:
            accession_variants.add(parts[1])
        if len(parts) >= 3:
            accession_variants.add(parts[2])
        if accession_variants & accessions:
            result.append(entry)
    return result


def filter_by_keyword(entries: List[oms.FASTAEntry], keyword: str) -> List[oms.FASTAEntry]:
    """Filter FASTA entries whose description or identifier contains the keyword."""
    keyword_lower = keyword.lower()
    result = []
    for entry in entries:
        if keyword_lower in entry.identifier.lower() or keyword_lower in entry.description.lower():
            result.append(entry)
    return result


def filter_by_length(
    entries: List[oms.FASTAEntry],
    min_length: Optional[int] = None,
    max_length: Optional[int] = None,
) -> List[oms.FASTAEntry]:
    """Filter FASTA entries by sequence length range."""
    result = []
    for entry in entries:
        seq_len = len(entry.sequence)
        if min_length is not None and seq_len < min_length:
            continue
        if max_length is not None and seq_len > max_length:
            continue
        result.append(entry)
    return result


def extract_subset(
    input_path: str,
    output_path: str,
    accessions_file: Optional[str] = None,
    keyword: Optional[str] = None,
    min_length: Optional[int] = None,
    max_length: Optional[int] = None,
) -> dict:
    """Extract a subset of proteins from a FASTA file based on the given criteria.

    Returns a dict with statistics about the extraction.
    """
    entries = load_fasta(input_path)
    total = len(entries)
    filtered = entries

    if accessions_file:
        with open(accessions_file) as fh:
            accessions = {line.strip() for line in fh if line.strip()}
        filtered = filter_by_accessions(filtered, accessions)

    if keyword:
        filtered = filter_by_keyword(filtered, keyword)

    if min_length is not None or max_length is not None:
        filtered = filter_by_length(filtered, min_length, max_length)

    save_fasta(filtered, output_path)
    return {"total_input": total, "total_output": len(filtered)}


@click.command(help="Extract proteins from a FASTA database by accession list, keyword, or length range.")
@click.option("--input", "input", required=True, help="Input FASTA file")
@click.option("--accessions", default=None, help="Text file with one accession per line")
@click.option("--keyword", default=None, help="Keyword to match in header/description")
@click.option("--min-length", type=int, default=None, help="Minimum sequence length")
@click.option("--max-length", type=int, default=None, help="Maximum sequence length")
@click.option("--output", required=True, help="Output FASTA file")
def main(input, accessions, keyword, min_length, max_length, output) -> None:
    if not accessions and not keyword and min_length is None and max_length is None:
        raise click.UsageError("At least one filter (--accessions, --keyword, --min-length, --max-length) is required.")

    stats = extract_subset(
        input, output, accessions, keyword, min_length, max_length
    )
    print(f"Extracted {stats['total_output']} / {stats['total_input']} proteins to {output}")


if __name__ == "__main__":
    main()
