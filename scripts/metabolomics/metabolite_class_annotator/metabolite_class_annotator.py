"""
Metabolite Class Annotator
============================
Annotate features with putative compound classes based on mass defect
analysis and elemental ratio heuristics.

Mass defect (fractional part of mass) and Kendrick mass defect are
used to classify features into broad compound families such as lipids,
peptides, sugars, and polyketides.

Usage
-----
    python metabolite_class_annotator.py --input features.tsv --output annotated.tsv
"""

import argparse
import csv
import sys

try:
    import pyopenms as oms  # noqa: F401
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

PROTON = 1.007276

# Compound class rules based on mass defect ranges (fractional mass)
# These are simplified heuristic boundaries.
CLASS_RULES = [
    {
        "name": "Lipid",
        "mass_defect_range": (0.0, 0.4),
        "mass_range": (200, 1200),
        "description": "Fatty acids, glycerolipids, sphingolipids",
    },
    {
        "name": "Peptide",
        "mass_defect_range": (0.0, 0.7),
        "mass_range": (400, 5000),
        "description": "Di- to oligopeptides",
    },
    {
        "name": "Sugar/Carbohydrate",
        "mass_defect_range": (0.0, 0.15),
        "mass_range": (100, 2000),
        "description": "Mono- to oligosaccharides",
    },
    {
        "name": "Polyketide",
        "mass_defect_range": (0.0, 0.35),
        "mass_range": (150, 800),
        "description": "Polyketide-derived natural products",
    },
    {
        "name": "Terpenoid",
        "mass_defect_range": (0.1, 0.5),
        "mass_range": (100, 800),
        "description": "Mono-, sesqui-, diterpenoids",
    },
    {
        "name": "Nucleoside",
        "mass_defect_range": (0.0, 0.2),
        "mass_range": (200, 600),
        "description": "Nucleosides and nucleotides",
    },
]


def compute_mass_defect(mass: float) -> float:
    """Compute the fractional mass defect.

    Parameters
    ----------
    mass:
        Monoisotopic mass in Da.

    Returns
    -------
    float
        Fractional part of the mass.
    """
    return mass - int(mass)


def compute_kendrick_mass_defect(mass: float, base_unit: float = 14.01565) -> float:
    """Compute the Kendrick mass defect (CH2-based).

    Parameters
    ----------
    mass:
        Monoisotopic mass in Da.
    base_unit:
        Kendrick base mass (default: CH2 = 14.01565 Da).

    Returns
    -------
    float
        Kendrick mass defect.
    """
    kendrick_mass = mass * (14.0 / base_unit)
    return round(kendrick_mass) - kendrick_mass


def annotate_class(mass: float) -> list[str]:
    """Annotate a mass with candidate compound classes.

    Parameters
    ----------
    mass:
        Neutral monoisotopic mass in Da.

    Returns
    -------
    list[str]
        List of candidate class names.
    """
    md = compute_mass_defect(mass)
    candidates = []

    for rule in CLASS_RULES:
        md_lo, md_hi = rule["mass_defect_range"]
        m_lo, m_hi = rule["mass_range"]
        if md_lo <= md <= md_hi and m_lo <= mass <= m_hi:
            candidates.append(rule["name"])

    return candidates if candidates else ["Unknown"]


def annotate_features(features: list[dict]) -> list[dict]:
    """Annotate a list of features with compound classes.

    Parameters
    ----------
    features:
        List of dicts with at least key ``mz``.

    Returns
    -------
    list[dict]
        Each feature augmented with mass_defect, kendrick_md, and compound_class.
    """
    results = []
    for feat in features:
        mz = float(feat["mz"])
        neutral_mass = mz - PROTON  # assume [M+H]+

        md = compute_mass_defect(neutral_mass)
        kmd = compute_kendrick_mass_defect(neutral_mass)
        classes = annotate_class(neutral_mass)

        feat_copy = dict(feat)
        feat_copy["neutral_mass"] = round(neutral_mass, 6)
        feat_copy["mass_defect"] = round(md, 6)
        feat_copy["kendrick_md"] = round(kmd, 6)
        feat_copy["compound_class"] = ";".join(classes)
        results.append(feat_copy)

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Annotate features with compound classes by mass defect analysis."
    )
    parser.add_argument("--input", required=True, metavar="FILE", help="Features TSV (must have mz column)")
    parser.add_argument("--output", required=True, metavar="FILE", help="Output annotated TSV")
    args = parser.parse_args()

    features = []
    with open(args.input) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            features.append(row)

    annotated = annotate_features(features)

    fieldnames = list(features[0].keys()) if features else ["mz"]
    fieldnames += ["neutral_mass", "mass_defect", "kendrick_md", "compound_class"]
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(annotated)

    print(f"Annotated {len(annotated)} features, written to {args.output}")


if __name__ == "__main__":
    main()
