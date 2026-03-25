"""
Mass Defect Filter
==================
Compute mass defect and Kendrick mass defect for features, then filter by MDF range.

Features:
- Compute mass defect (fractional part of exact mass)
- Kendrick mass defect with configurable base (e.g. CH2)
- Filter features by mass defect range
- TSV input/output

Usage
-----
    python mass_defect_filter.py --input features.tsv --mdf-min 0.1 --mdf-max 0.3
    python mass_defect_filter.py --input features.tsv --mdf-min 0.1 --mdf-max 0.3 --kendrick-base CH2 \
        --output filtered.tsv
"""

import argparse
import csv
import sys

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def compute_mass_defect(exact_mass: float) -> float:
    """Compute mass defect (fractional part of exact mass).

    Parameters
    ----------
    exact_mass : float
        Exact monoisotopic mass.

    Returns
    -------
    float
        Mass defect value.
    """
    return exact_mass - int(exact_mass)


def compute_kendrick_mass(exact_mass: float, kendrick_base: str = "CH2") -> float:
    """Compute Kendrick mass using given base.

    Parameters
    ----------
    exact_mass : float
        Exact monoisotopic mass.
    kendrick_base : str
        Molecular formula for the Kendrick base (default "CH2").

    Returns
    -------
    float
        Kendrick mass.
    """
    formula = oms.EmpiricalFormula(kendrick_base)
    exact_base = formula.getMonoWeight()
    nominal_base = round(exact_base)
    if exact_base == 0:
        return exact_mass
    return exact_mass * (nominal_base / exact_base)


def compute_kendrick_mass_defect(exact_mass: float, kendrick_base: str = "CH2") -> float:
    """Compute Kendrick mass defect.

    Parameters
    ----------
    exact_mass : float
        Exact monoisotopic mass.
    kendrick_base : str
        Molecular formula for the Kendrick base.

    Returns
    -------
    float
        Kendrick mass defect.
    """
    km = compute_kendrick_mass(exact_mass, kendrick_base)
    return round(km) - km


def filter_by_mass_defect(
    features: list[dict],
    mdf_min: float = 0.0,
    mdf_max: float = 1.0,
    kendrick_base: str = "CH2",
) -> list[dict]:
    """Filter features by mass defect range and compute Kendrick mass defect.

    Parameters
    ----------
    features : list[dict]
        List of feature dicts, each must have an 'exact_mass' key.
    mdf_min : float
        Minimum mass defect (inclusive).
    mdf_max : float
        Maximum mass defect (inclusive).
    kendrick_base : str
        Molecular formula for Kendrick base.

    Returns
    -------
    list[dict]
        Filtered features with additional mass_defect and kendrick_mass_defect keys.
    """
    results = []
    for feat in features:
        mass = float(feat["exact_mass"])
        md = compute_mass_defect(mass)
        kmd = compute_kendrick_mass_defect(mass, kendrick_base)

        if mdf_min <= md <= mdf_max:
            result = dict(feat)
            result["mass_defect"] = round(md, 6)
            result["kendrick_mass_defect"] = round(kmd, 6)
            results.append(result)

    return results


def read_features_tsv(input_path: str) -> list[dict]:
    """Read features from a TSV file.

    Parameters
    ----------
    input_path : str
        Path to TSV file with at least an 'exact_mass' column.

    Returns
    -------
    list[dict]
        List of feature dictionaries.
    """
    features = []
    with open(input_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            features.append(row)
    return features


def write_tsv(results: list[dict], output_path: str) -> None:
    """Write filtered features to TSV file.

    Parameters
    ----------
    results : list[dict]
        List of result dictionaries.
    output_path : str
        Path to output TSV file.
    """
    if not results:
        with open(output_path, "w") as f:
            f.write("exact_mass\tmass_defect\tkendrick_mass_defect\n")
        return

    fieldnames = list(results[0].keys())
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


def main():
    parser = argparse.ArgumentParser(
        description="Compute mass defect and Kendrick mass defect, filter features."
    )
    parser.add_argument("--input", required=True, help="Input TSV with exact_mass column")
    parser.add_argument("--mdf-min", type=float, default=0.0, help="Minimum mass defect (default: 0.0)")
    parser.add_argument("--mdf-max", type=float, default=1.0, help="Maximum mass defect (default: 1.0)")
    parser.add_argument("--kendrick-base", default="CH2", help="Kendrick base formula (default: CH2)")
    parser.add_argument("--output", default=None, help="Output TSV file path (default: print to stdout)")
    args = parser.parse_args()

    features = read_features_tsv(args.input)
    results = filter_by_mass_defect(features, args.mdf_min, args.mdf_max, args.kendrick_base)

    if args.output:
        write_tsv(results, args.output)
        print(f"Wrote {len(results)} filtered features to {args.output}")
    else:
        if results:
            print("\t".join(results[0].keys()))
            for r in results:
                print("\t".join(str(v) for v in r.values()))
        else:
            print("No features matched the mass defect filter.")


if __name__ == "__main__":
    main()
