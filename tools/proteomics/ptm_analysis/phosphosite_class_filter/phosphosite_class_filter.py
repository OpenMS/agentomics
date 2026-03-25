"""
Phosphosite Class Filter
========================
Classify phosphosites into Class I/II/III by localization probability.
Report enrichment efficiency.

Class I:   localization_prob >= class1_threshold (default 0.75)
Class II:  0.50 <= localization_prob < class1_threshold
Class III: localization_prob < 0.50

Usage
-----
    python phosphosite_class_filter.py --input phosphosites.tsv --class1-threshold 0.75 --output classified.tsv
"""

import csv
import sys
from typing import Dict, List, Tuple

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


CLASS1_DEFAULT_THRESHOLD = 0.75
CLASS2_LOWER = 0.50


def classify_phosphosite(localization_prob: float, class1_threshold: float = CLASS1_DEFAULT_THRESHOLD) -> str:
    """Classify a phosphosite into Class I, II, or III based on localization probability.

    Parameters
    ----------
    localization_prob:
        Localization probability (0.0 to 1.0).
    class1_threshold:
        Minimum probability for Class I classification.

    Returns
    -------
    str
        One of 'Class I', 'Class II', 'Class III'.
    """
    if localization_prob >= class1_threshold:
        return "Class I"
    elif localization_prob >= CLASS2_LOWER:
        return "Class II"
    else:
        return "Class III"


def validate_modification(modification: str) -> bool:
    """Check if a modification string refers to a phosphorylation using ModificationsDB.

    Parameters
    ----------
    modification:
        Modification name string.

    Returns
    -------
    bool
        True if the modification is recognized as phosphorylation-related.
    """
    phospho_keywords = ["phospho", "phos", "79.966"]
    mod_lower = modification.lower()
    for kw in phospho_keywords:
        if kw in mod_lower:
            return True
    # Try to look up via ModificationsDB
    try:
        mod_db = oms.ModificationsDB()
        mod_obj = oms.ResidueModification()
        mod_db.getModification(modification, mod_obj)
        diff = mod_obj.getDiffMonoMass()
        # Phosphorylation adds ~79.966 Da
        if abs(diff - 79.966) < 0.1:
            return True
    except Exception:
        pass
    return False


def validate_phosphopeptide(sequence: str) -> bool:
    """Validate that a sequence can be parsed as an AASequence.

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


def classify_phosphosites(
    rows: List[Dict[str, str]],
    class1_threshold: float = CLASS1_DEFAULT_THRESHOLD,
) -> Tuple[List[Dict[str, str]], Dict[str, int]]:
    """Classify a list of phosphosite rows and compute summary statistics.

    Parameters
    ----------
    rows:
        List of dicts with keys: peptide, protein, site, localization_prob, modification.
    class1_threshold:
        Minimum probability for Class I classification.

    Returns
    -------
    tuple
        (classified_rows, summary) where summary has counts per class.
    """
    summary: Dict[str, int] = {"Class I": 0, "Class II": 0, "Class III": 0, "total": 0}
    classified = []

    for row in rows:
        prob = float(row["localization_prob"])
        site_class = classify_phosphosite(prob, class1_threshold)
        new_row = dict(row)
        new_row["site_class"] = site_class
        new_row["valid_peptide"] = str(validate_phosphopeptide(row["peptide"]))
        classified.append(new_row)
        summary[site_class] += 1
        summary["total"] += 1

    return classified, summary


def compute_enrichment_efficiency(summary: Dict[str, int]) -> float:
    """Compute enrichment efficiency as fraction of Class I sites.

    Parameters
    ----------
    summary:
        Dictionary with class counts and total.

    Returns
    -------
    float
        Fraction of Class I sites (0.0 to 1.0).
    """
    if summary["total"] == 0:
        return 0.0
    return summary["Class I"] / summary["total"]


def read_input(input_path: str) -> List[Dict[str, str]]:
    """Read phosphosite TSV input file.

    Parameters
    ----------
    input_path:
        Path to input TSV file.

    Returns
    -------
    list
        List of row dictionaries.
    """
    rows = []
    with open(input_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def write_output(output_path: str, classified_rows: List[Dict[str, str]], summary: Dict[str, int]) -> None:
    """Write classified phosphosites and summary to TSV.

    Parameters
    ----------
    output_path:
        Path to output TSV file.
    classified_rows:
        List of classified row dictionaries.
    summary:
        Summary statistics dictionary.
    """
    if not classified_rows:
        return

    fieldnames = list(classified_rows[0].keys())
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(classified_rows)

    # Write summary to a companion file
    summary_path = output_path.replace(".tsv", "_summary.tsv")
    enrichment = compute_enrichment_efficiency(summary)
    with open(summary_path, "w", newline="") as f:
        f.write("class\tcount\n")
        for cls in ["Class I", "Class II", "Class III"]:
            f.write(f"{cls}\t{summary[cls]}\n")
        f.write(f"total\t{summary['total']}\n")
        f.write(f"enrichment_efficiency\t{enrichment:.4f}\n")


@click.command(help="Classify phosphosites into Class I/II/III by localization probability.")
@click.option("--input", "input", required=True, help="Input phosphosites TSV file")
@click.option(
    "--class1-threshold", type=float, default=CLASS1_DEFAULT_THRESHOLD,
    help=f"Minimum localization probability for Class I (default: {CLASS1_DEFAULT_THRESHOLD})",
)
@click.option("--output", required=True, help="Output classified TSV file")
def main(input, class1_threshold, output):
    rows = read_input(input)
    classified, summary = classify_phosphosites(rows, class1_threshold)
    write_output(output, classified, summary)

    enrichment = compute_enrichment_efficiency(summary)
    print(f"Total phosphosites: {summary['total']}")
    print(f"  Class I:  {summary['Class I']}")
    print(f"  Class II: {summary['Class II']}")
    print(f"  Class III: {summary['Class III']}")
    print(f"Enrichment efficiency (Class I fraction): {enrichment:.4f}")


if __name__ == "__main__":
    main()
