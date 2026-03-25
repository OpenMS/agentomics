"""
In-Source Fragmentation Detector
=================================
Detect in-source fragmentation (ISF) artifacts by identifying coeluting features
whose mass difference matches common neutral losses (H2O, CO2, NH3, etc.).

Usage
-----
    python isf_detector.py --input features.tsv --rt-tolerance 3 --output isf_annotated.tsv
"""

import csv
import math
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


# Common neutral losses and their exact masses
NEUTRAL_LOSSES = {
    "H2O": oms.EmpiricalFormula("H2O").getMonoWeight(),
    "CO2": oms.EmpiricalFormula("CO2").getMonoWeight(),
    "NH3": oms.EmpiricalFormula("NH3").getMonoWeight(),
    "CO": oms.EmpiricalFormula("CO").getMonoWeight(),
    "HCOOH": oms.EmpiricalFormula("CH2O2").getMonoWeight(),
    "CH3OH": oms.EmpiricalFormula("CH4O").getMonoWeight(),
}


def get_neutral_loss_masses() -> dict:
    """Return the neutral loss name-to-mass mapping.

    Returns
    -------
    dict mapping loss names to exact masses.
    """
    return dict(NEUTRAL_LOSSES)


def pearson_correlation(x: list, y: list) -> float:
    """Compute Pearson correlation coefficient between two lists.

    Parameters
    ----------
    x, y:
        Lists of numeric values of equal length.

    Returns
    -------
    float: Pearson r value, or 0.0 if computation is not possible.
    """
    n = len(x)
    if n < 2:
        return 0.0

    mean_x = sum(x) / n
    mean_y = sum(y) / n

    cov = sum((xi - mean_x) * (yi - mean_y) for xi, yi in zip(x, y))
    std_x = math.sqrt(sum((xi - mean_x) ** 2 for xi in x))
    std_y = math.sqrt(sum((yi - mean_y) ** 2 for yi in y))

    if std_x == 0 or std_y == 0:
        return 0.0

    return cov / (std_x * std_y)


def detect_isf_pairs(
    features: list,
    rt_tolerance: float = 3.0,
    mass_tolerance_da: float = 0.01,
) -> list:
    """Detect potential in-source fragmentation pairs.

    Parameters
    ----------
    features:
        List of dicts with keys: id, mz, rt, intensity.
    rt_tolerance:
        Maximum RT difference (seconds) for coelution.
    mass_tolerance_da:
        Mass tolerance in Daltons for neutral loss matching.

    Returns
    -------
    list of dicts describing ISF pair annotations.
    """
    pairs = []
    n = len(features)

    for i in range(n):
        for j in range(i + 1, n):
            fi = features[i]
            fj = features[j]

            rt_diff = abs(float(fi["rt"]) - float(fj["rt"]))
            if rt_diff > rt_tolerance:
                continue

            mass_diff = abs(float(fi["mz"]) - float(fj["mz"]))

            for loss_name, loss_mass in NEUTRAL_LOSSES.items():
                if abs(mass_diff - loss_mass) <= mass_tolerance_da:
                    # Heavier feature is the precursor
                    if float(fi["mz"]) > float(fj["mz"]):
                        precursor, fragment = fi, fj
                    else:
                        precursor, fragment = fj, fi

                    pairs.append({
                        "precursor_id": precursor["id"],
                        "precursor_mz": float(precursor["mz"]),
                        "fragment_id": fragment["id"],
                        "fragment_mz": float(fragment["mz"]),
                        "neutral_loss": loss_name,
                        "mass_diff": round(mass_diff, 6),
                        "rt_diff": round(rt_diff, 2),
                        "isf_candidate": True,
                    })
                    break  # Only assign one neutral loss per pair

    return pairs


def annotate_features(features: list, isf_pairs: list) -> list:
    """Add ISF annotation columns to the feature list.

    Parameters
    ----------
    features:
        Original feature list.
    isf_pairs:
        ISF pair annotations from detect_isf_pairs.

    Returns
    -------
    list of dicts with added isf_flag, isf_role, isf_partner, isf_loss columns.
    """
    # Build lookup for partner info
    partner_info = {}
    for p in isf_pairs:
        partner_info[p["fragment_id"]] = {
            "role": "fragment",
            "partner": p["precursor_id"],
            "loss": p["neutral_loss"],
        }
        if p["precursor_id"] not in partner_info:
            partner_info[p["precursor_id"]] = {
                "role": "precursor",
                "partner": p["fragment_id"],
                "loss": p["neutral_loss"],
            }

    annotated = []
    for f in features:
        row = dict(f)
        fid = f["id"]
        if fid in partner_info:
            row["isf_flag"] = True
            row["isf_role"] = partner_info[fid]["role"]
            row["isf_partner"] = partner_info[fid]["partner"]
            row["isf_loss"] = partner_info[fid]["loss"]
        else:
            row["isf_flag"] = False
            row["isf_role"] = ""
            row["isf_partner"] = ""
            row["isf_loss"] = ""
        annotated.append(row)

    return annotated


@click.command()
@click.option("--input", "input_file", required=True, help="TSV with columns: id, mz, rt, intensity.")
@click.option("--rt-tolerance", type=float, default=3.0,
              help="RT tolerance in seconds (default: 3).")
@click.option("--mass-tolerance", type=float, default=0.01,
              help="Mass tolerance in Da for neutral loss matching (default: 0.01).")
@click.option("--output", required=True, help="Output TSV with ISF annotations.")
def main(input_file, rt_tolerance, mass_tolerance, output) -> None:
    """CLI entry point."""
    features = []
    with open(input_file, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            features.append(row)

    isf_pairs = detect_isf_pairs(features, rt_tolerance=rt_tolerance,
                                 mass_tolerance_da=mass_tolerance)
    annotated = annotate_features(features, isf_pairs)

    fieldnames = list(annotated[0].keys()) if annotated else []
    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(annotated)

    n_isf = sum(1 for a in annotated if a["isf_flag"])
    print(f"Wrote {len(annotated)} features ({n_isf} ISF-flagged) to {output}")


if __name__ == "__main__":
    main()
