"""
FASTA Decoy Validator
=====================
Check if a FASTA database contains decoy sequences and validate prefix consistency.

Usage
-----
    python fasta_decoy_validator.py --input db.fasta --decoy-prefix DECOY_ --output validation.json
"""

import argparse
import json
import sys
from typing import List

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


def validate_decoys(
    input_path: str,
    decoy_prefix: str = "DECOY_",
) -> dict:
    """Validate decoy sequences in a FASTA database.

    Returns a dict with validation results including:
    - total_entries: total number of proteins
    - target_count: number of target entries
    - decoy_count: number of decoy entries
    - has_decoys: whether decoys were found
    - decoy_ratio: ratio of decoys to targets
    - prefix_consistent: whether all decoys use the expected prefix
    - alternative_prefixes: any other prefixes detected
    - reversed_match_count: number of decoys whose sequence is the reverse of a target
    """
    entries = load_fasta(input_path)
    total = len(entries)

    common_prefixes = ["DECOY_", "REV_", "rev_", "decoy_", "REVERSED_", "XXX_"]
    if decoy_prefix not in common_prefixes:
        common_prefixes.insert(0, decoy_prefix)

    target_entries = []
    decoy_entries = []
    prefix_counts = {}

    for entry in entries:
        identifier = entry.identifier
        is_decoy = False
        for prefix in common_prefixes:
            if identifier.startswith(prefix):
                prefix_counts[prefix] = prefix_counts.get(prefix, 0) + 1
                is_decoy = True
                break
        if is_decoy:
            decoy_entries.append(entry)
        else:
            target_entries.append(entry)

    target_count = len(target_entries)
    decoy_count = len(decoy_entries)
    expected_prefix_count = prefix_counts.get(decoy_prefix, 0)

    # Check for prefix consistency
    alternative_prefixes = {p: c for p, c in prefix_counts.items() if p != decoy_prefix and c > 0}
    prefix_consistent = decoy_count == expected_prefix_count

    # Check if decoys are reversed versions of targets
    target_seqs = {e.sequence for e in target_entries}
    reversed_match = 0
    for entry in decoy_entries:
        if entry.sequence[::-1] in target_seqs:
            reversed_match += 1

    decoy_ratio = decoy_count / target_count if target_count > 0 else 0.0

    return {
        "total_entries": total,
        "target_count": target_count,
        "decoy_count": decoy_count,
        "has_decoys": decoy_count > 0,
        "decoy_ratio": round(decoy_ratio, 4),
        "expected_prefix": decoy_prefix,
        "expected_prefix_count": expected_prefix_count,
        "prefix_consistent": prefix_consistent,
        "alternative_prefixes": alternative_prefixes,
        "reversed_match_count": reversed_match,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate decoy sequences in a FASTA database."
    )
    parser.add_argument("--input", required=True, help="Input FASTA file")
    parser.add_argument("--decoy-prefix", default="DECOY_", help="Expected decoy prefix (default: DECOY_)")
    parser.add_argument("--output", default=None, help="Output JSON file (default: stdout)")
    args = parser.parse_args()

    result = validate_decoys(args.input, args.decoy_prefix)
    output = json.dumps(result, indent=2)

    if args.output:
        with open(args.output, "w") as fh:
            fh.write(output + "\n")
        print(f"Validation results written to {args.output}")
    else:
        print(output)


if __name__ == "__main__":
    main()
