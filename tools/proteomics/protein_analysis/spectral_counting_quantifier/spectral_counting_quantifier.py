"""
Spectral Counting Quantifier
==============================
Calculate protein abundances from spectral counts using emPAI and NSAF methods.

Features
--------
- emPAI (exponentially modified Protein Abundance Index)
- NSAF (Normalized Spectral Abundance Factor)
- Read peptide-spectrum counts from TSV input
- In-silico digestion for observable peptide count (emPAI)

Usage
-----
    python spectral_counting_quantifier.py --input counts.tsv --fasta db.fasta --method nsaf --output out.tsv
    python spectral_counting_quantifier.py --input counts.tsv --fasta db.fasta --method empai --output out.tsv
"""

import argparse
import csv
import json
import sys

try:
    import pyopenms as oms  # noqa: F401
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def load_fasta_proteins(fasta_path: str) -> dict:
    """Load protein sequences from FASTA.

    Parameters
    ----------
    fasta_path : str
        Path to FASTA file.

    Returns
    -------
    dict
        Mapping of accession to sequence.
    """
    entries = []
    oms.FASTAFile().load(fasta_path, entries)
    return {e.identifier: e.sequence for e in entries}


def count_observable_peptides(sequence: str, enzyme: str = "Trypsin",
                              min_length: int = 6, max_length: int = 40) -> int:
    """Count observable tryptic peptides for a protein sequence.

    Parameters
    ----------
    sequence : str
        Protein sequence.
    enzyme : str
        Enzyme name.
    min_length : int
        Minimum peptide length.
    max_length : int
        Maximum peptide length.

    Returns
    -------
    int
        Number of observable peptides.
    """
    aa_seq = oms.AASequence.fromString(sequence)
    digest = oms.ProteaseDigestion()
    digest.setEnzyme(enzyme)
    digest.setMissedCleavages(0)
    peptides = []
    digest.digest(aa_seq, peptides, min_length, max_length)
    return len(peptides)


def calculate_empai(protein_data: dict, proteins: dict, enzyme: str = "Trypsin") -> list:
    """Calculate emPAI for each protein.

    Parameters
    ----------
    protein_data : dict
        Mapping of protein accession to dict with 'spectral_count' and 'observed_peptides'.
    proteins : dict
        Mapping of accession to sequence.
    enzyme : str
        Enzyme name.

    Returns
    -------
    list
        List of dicts with emPAI values.
    """
    results = []
    for accession, data in protein_data.items():
        seq = proteins.get(accession, "")
        n_observable = count_observable_peptides(seq, enzyme) if seq else 1
        n_observed = data.get("observed_peptides", 0)
        pai = n_observed / max(n_observable, 1)
        empai = 10 ** pai - 1

        results.append({
            "accession": accession,
            "spectral_count": data.get("spectral_count", 0),
            "observed_peptides": n_observed,
            "observable_peptides": n_observable,
            "pai": round(pai, 6),
            "empai": round(empai, 6),
        })
    return results


def calculate_nsaf(protein_data: dict, proteins: dict) -> list:
    """Calculate NSAF for each protein.

    Parameters
    ----------
    protein_data : dict
        Mapping of protein accession to dict with 'spectral_count'.
    proteins : dict
        Mapping of accession to sequence.

    Returns
    -------
    list
        List of dicts with NSAF values.
    """
    # Calculate SAF = SpC / Length for each protein
    saf_values = {}
    for accession, data in protein_data.items():
        seq = proteins.get(accession, "")
        length = len(seq) if seq else 1
        spc = data.get("spectral_count", 0)
        saf_values[accession] = spc / max(length, 1)

    total_saf = sum(saf_values.values())

    results = []
    for accession, data in protein_data.items():
        seq = proteins.get(accession, "")
        length = len(seq) if seq else 0
        saf = saf_values[accession]
        nsaf = saf / total_saf if total_saf > 0 else 0.0

        results.append({
            "accession": accession,
            "spectral_count": data.get("spectral_count", 0),
            "protein_length": length,
            "saf": round(saf, 8),
            "nsaf": round(nsaf, 8),
        })
    return results


def load_peptide_counts(input_path: str) -> dict:
    """Load peptide spectral counts from TSV and aggregate per protein.

    Parameters
    ----------
    input_path : str
        Path to TSV file with columns: protein, peptide, spectral_count.

    Returns
    -------
    dict
        Mapping of protein accession to aggregated counts.
    """
    protein_data = {}
    with open(input_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            protein = row.get("protein", "").strip()
            spc = int(row.get("spectral_count", 0))
            if protein not in protein_data:
                protein_data[protein] = {"spectral_count": 0, "observed_peptides": 0, "peptides": set()}
            protein_data[protein]["spectral_count"] += spc
            peptide = row.get("peptide", "").strip()
            if peptide and peptide not in protein_data[protein]["peptides"]:
                protein_data[protein]["peptides"].add(peptide)
                protein_data[protein]["observed_peptides"] += 1

    # Remove set (not serializable)
    for v in protein_data.values():
        del v["peptides"]
    return protein_data


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(description="Calculate protein abundances from spectral counts.")
    parser.add_argument("--input", required=True, help="TSV with protein, peptide, spectral_count columns.")
    parser.add_argument("--fasta", required=True, help="Protein FASTA database.")
    parser.add_argument("--method", choices=["empai", "nsaf"], default="nsaf", help="Quantification method.")
    parser.add_argument("--enzyme", default="Trypsin", help="Enzyme for emPAI (default: Trypsin).")
    parser.add_argument("--output", help="Output file (.tsv or .json).")
    args = parser.parse_args()

    proteins = load_fasta_proteins(args.fasta)
    protein_data = load_peptide_counts(args.input)

    if args.method == "empai":
        results = calculate_empai(protein_data, proteins, args.enzyme)
    else:
        results = calculate_nsaf(protein_data, proteins)

    if args.output:
        if args.output.endswith(".json"):
            with open(args.output, "w") as fh:
                json.dump(results, fh, indent=2)
        else:
            with open(args.output, "w", newline="") as fh:
                fieldnames = list(results[0].keys()) if results else []
                writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
                writer.writeheader()
                writer.writerows(results)
        print(f"Results written to {args.output}")
    else:
        for r in results:
            print(json.dumps(r))


if __name__ == "__main__":
    main()
