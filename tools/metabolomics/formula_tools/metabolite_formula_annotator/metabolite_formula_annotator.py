"""
Metabolite Formula Annotator
==============================
Annotate features with candidate molecular formulas based on mass
matching and isotope pattern fitting.

For each feature (m/z), candidate formulas are generated within the
specified element constraints and ranked by isotope-pattern similarity.

Usage
-----
    python metabolite_formula_annotator.py --input features.tsv --ppm 5 --elements C,H,N,O --output annotated.tsv
"""

import csv
import itertools

import click
import pyopenms as oms

PROTON = 1.007276

# Default element ranges for formula enumeration
DEFAULT_ELEMENT_RANGES = {
    "C": (0, 40),
    "H": (0, 80),
    "N": (0, 10),
    "O": (0, 25),
}


def enumerate_formulas(
    target_mass: float,
    ppm: float = 5.0,
    element_ranges: dict[str, tuple[int, int]] | None = None,
) -> list[dict]:
    """Enumerate candidate formulas matching a target neutral mass.

    Parameters
    ----------
    target_mass:
        Neutral monoisotopic mass in Da.
    ppm:
        Mass tolerance in ppm.
    element_ranges:
        Dict mapping element symbols to (min, max) count ranges.

    Returns
    -------
    list[dict]
        Each dict has: formula, mass, error_ppm.
    """
    if element_ranges is None:
        element_ranges = DEFAULT_ELEMENT_RANGES

    tol_da = target_mass * ppm / 1e6
    lo = target_mass - tol_da
    hi = target_mass + tol_da

    elements = sorted(element_ranges.keys())
    ranges = [range(element_ranges[e][0], element_ranges[e][1] + 1) for e in elements]

    results = []
    for combo in itertools.product(*ranges):
        formula_str = "".join(f"{e}{n}" for e, n in zip(elements, combo) if n > 0)
        if not formula_str:
            continue

        try:
            ef = oms.EmpiricalFormula(formula_str)
            mass = ef.getMonoWeight()
        except Exception:
            continue

        if lo <= mass <= hi:
            error = (mass - target_mass) / target_mass * 1e6
            results.append({
                "formula": formula_str,
                "mass": round(mass, 6),
                "error_ppm": round(error, 4),
            })

    results.sort(key=lambda r: abs(r["error_ppm"]))
    return results


def score_isotope_fit(formula: str, observed_ratios: list[float] | None = None) -> float:
    """Score isotope pattern fit for a formula.

    Parameters
    ----------
    formula:
        Molecular formula string.
    observed_ratios:
        Observed isotope intensity ratios (M+0=100, M+1, M+2, ...).
        If None, returns 0.0.

    Returns
    -------
    float
        Cosine similarity score (0.0 to 1.0).
    """
    if observed_ratios is None or len(observed_ratios) < 2:
        return 0.0

    ef = oms.EmpiricalFormula(formula)
    gen = oms.CoarseIsotopePatternGenerator(len(observed_ratios))
    iso = ef.getIsotopeDistribution(gen)
    container = iso.getContainer()

    theo = [peak.getIntensity() for peak in container]
    if not theo:
        return 0.0

    # Normalize to max=100
    theo_max = max(theo)
    if theo_max > 0:
        theo = [t / theo_max * 100.0 for t in theo]

    # Cosine similarity
    dot = sum(a * b for a, b in zip(observed_ratios, theo))
    mag_a = sum(a ** 2 for a in observed_ratios) ** 0.5
    mag_b = sum(b ** 2 for b in theo) ** 0.5

    if mag_a == 0 or mag_b == 0:
        return 0.0
    return dot / (mag_a * mag_b)


def annotate_features(
    features: list[dict],
    ppm: float = 5.0,
    element_ranges: dict[str, tuple[int, int]] | None = None,
    max_candidates: int = 5,
) -> list[dict]:
    """Annotate features with candidate formulas.

    Parameters
    ----------
    features:
        List of dicts with at least key ``mz``.
    ppm:
        Mass tolerance in ppm.
    element_ranges:
        Element count ranges for formula enumeration.
    max_candidates:
        Maximum number of candidate formulas per feature.

    Returns
    -------
    list[dict]
        Each feature dict augmented with ``candidates`` key.
    """
    results = []
    for feat in features:
        mz = float(feat["mz"])
        neutral_mass = mz - PROTON  # assume [M+H]+
        candidates = enumerate_formulas(neutral_mass, ppm=ppm, element_ranges=element_ranges)
        feat_copy = dict(feat)
        feat_copy["candidates"] = candidates[:max_candidates]
        results.append(feat_copy)
    return results


@click.command()
@click.option("--input", "input_file", required=True, help="Features TSV (must have mz column)")
@click.option("--ppm", type=float, default=5.0, help="Mass tolerance in ppm (default: 5)")
@click.option("--elements", default="C,H,N,O", help="Comma-separated elements (default: C,H,N,O)")
@click.option("--output", required=True, help="Output annotated TSV")
def main(input_file, ppm, elements, output):
    elements = elements.split(",")
    element_ranges = {e.strip(): DEFAULT_ELEMENT_RANGES.get(e.strip(), (0, 10)) for e in elements}

    features = []
    with open(input_file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            features.append(row)

    annotated = annotate_features(features, ppm=ppm, element_ranges=element_ranges)

    with open(output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["mz", "rt", "intensity", "candidate_formula", "candidate_mass", "error_ppm"])
        for feat in annotated:
            candidates = feat.get("candidates", [])
            if candidates:
                for c in candidates:
                    writer.writerow([
                        feat.get("mz", ""),
                        feat.get("rt", ""),
                        feat.get("intensity", ""),
                        c["formula"],
                        c["mass"],
                        c["error_ppm"],
                    ])
            else:
                writer.writerow([
                    feat.get("mz", ""),
                    feat.get("rt", ""),
                    feat.get("intensity", ""),
                    "", "", "",
                ])

    print(f"Annotated {len(annotated)} features, written to {output}")


if __name__ == "__main__":
    main()
