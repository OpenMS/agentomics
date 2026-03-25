"""
mzTab Summarizer
================
Parse an mzTab file and extract summary statistics.

Usage
-----
    python mztab_summarizer.py --input results.mzTab --output summary.tsv
"""

import csv
import sys
from collections import Counter

import click

try:
    import pyopenms as oms  # noqa: F401
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def parse_mztab(input_path: str) -> dict:
    """Parse an mzTab file manually and return its sections.

    Returns a dict with keys: metadata, protein_header, proteins,
    peptide_header, peptides, psm_header, psms.
    """
    metadata = {}
    protein_header = []
    proteins = []
    peptide_header = []
    peptides = []
    psm_header = []
    psms = []

    with open(input_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            prefix = parts[0]

            if prefix == "MTD":
                if len(parts) >= 3:
                    metadata[parts[1]] = parts[2]
            elif prefix == "PRH":
                protein_header = parts[1:]
            elif prefix == "PRT":
                proteins.append(parts[1:])
            elif prefix == "PEH":
                peptide_header = parts[1:]
            elif prefix == "PEP":
                peptides.append(parts[1:])
            elif prefix == "PSH":
                psm_header = parts[1:]
            elif prefix == "PSM":
                psms.append(parts[1:])

    return {
        "metadata": metadata,
        "protein_header": protein_header,
        "proteins": proteins,
        "peptide_header": peptide_header,
        "peptides": peptides,
        "psm_header": psm_header,
        "psms": psms,
    }


def summarize_mztab(input_path: str) -> dict:
    """Compute summary statistics from an mzTab file.

    Returns a dict with:
    - protein_count: total proteins
    - peptide_count: total peptides
    - psm_count: total PSMs
    - unique_proteins: unique protein accessions
    - unique_peptides: unique peptide sequences
    - modification_counts: modification frequency
    - search_engines: search engines used
    """
    data = parse_mztab(input_path)

    # Protein stats
    protein_accessions = set()
    if data["protein_header"]:
        acc_idx = data["protein_header"].index("accession") if "accession" in data["protein_header"] else 0
        for row in data["proteins"]:
            if len(row) > acc_idx:
                protein_accessions.add(row[acc_idx])

    # Peptide stats
    unique_peptide_seqs = set()
    mod_counter: Counter = Counter()
    if data["peptide_header"]:
        seq_idx = data["peptide_header"].index("sequence") if "sequence" in data["peptide_header"] else 0
        mod_idx = (
            data["peptide_header"].index("modifications")
            if "modifications" in data["peptide_header"]
            else None
        )
        for row in data["peptides"]:
            if len(row) > seq_idx:
                unique_peptide_seqs.add(row[seq_idx])
            if mod_idx is not None and len(row) > mod_idx and row[mod_idx] != "null":
                mod_counter[row[mod_idx]] += 1

    # PSM stats
    unique_psm_peptides = set()
    if data["psm_header"]:
        seq_idx = data["psm_header"].index("sequence") if "sequence" in data["psm_header"] else 0
        for row in data["psms"]:
            if len(row) > seq_idx:
                unique_psm_peptides.add(row[seq_idx])

    # Search engine info
    search_engines = [v for k, v in data["metadata"].items() if "search_engine" in k.lower()]

    return {
        "protein_count": len(data["proteins"]),
        "peptide_count": len(data["peptides"]),
        "psm_count": len(data["psms"]),
        "unique_protein_accessions": len(protein_accessions),
        "unique_peptide_sequences": len(unique_peptide_seqs),
        "unique_psm_peptides": len(unique_psm_peptides),
        "modification_counts": dict(mod_counter),
        "search_engines": search_engines,
        "metadata_entries": len(data["metadata"]),
    }


def write_summary_tsv(summary: dict, output_path: str) -> None:
    """Write summary statistics to a TSV file."""
    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key, value in summary.items():
            if isinstance(value, dict):
                for k2, v2 in value.items():
                    writer.writerow([f"{key}:{k2}", v2])
            elif isinstance(value, list):
                writer.writerow([key, ";".join(str(v) for v in value)])
            else:
                writer.writerow([key, value])


@click.command(help="Parse mzTab and extract summary statistics.")
@click.option("--input", "input", required=True, help="Input mzTab file")
@click.option("--output", default=None, help="Output summary TSV file (default: stdout)")
def main(input, output) -> None:
    summary = summarize_mztab(input)

    if output:
        write_summary_tsv(summary, output)
        print(f"Summary written to {output}")
    else:
        for key, value in summary.items():
            print(f"{key}: {value}")


if __name__ == "__main__":
    main()
