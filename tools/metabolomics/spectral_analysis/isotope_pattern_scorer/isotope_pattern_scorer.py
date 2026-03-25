"""
Isotope Pattern Scorer
=======================
Score observed isotope patterns against theoretical patterns computed
from a molecular formula using pyopenms.

The score is a cosine similarity between observed and theoretical
isotope intensity ratios.

Usage
-----
    python isotope_pattern_scorer.py --observed "180.063:100,181.067:6.5" --formula C6H12O6 --output fit.json
"""

import json
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def parse_observed(observed_str: str) -> list[tuple[float, float]]:
    """Parse observed isotope pattern from a string.

    Parameters
    ----------
    observed_str:
        Comma-separated ``mz:intensity`` pairs, e.g.
        ``"180.063:100,181.067:6.5"``.

    Returns
    -------
    list[tuple[float, float]]
        List of (mz, intensity) tuples.
    """
    peaks = []
    for pair in observed_str.split(","):
        pair = pair.strip()
        if ":" in pair:
            mz_str, int_str = pair.split(":", 1)
            peaks.append((float(mz_str), float(int_str)))
    return peaks


def get_theoretical_pattern(formula: str, n_peaks: int = 5) -> list[tuple[float, float]]:
    """Compute theoretical isotope pattern for a formula.

    Parameters
    ----------
    formula:
        Molecular formula string.
    n_peaks:
        Number of isotope peaks to generate.

    Returns
    -------
    list[tuple[float, float]]
        List of (mz, relative_intensity) tuples, normalized to max=100.
    """
    ef = oms.EmpiricalFormula(formula)
    gen = oms.CoarseIsotopePatternGenerator(n_peaks)
    iso = ef.getIsotopeDistribution(gen)
    container = iso.getContainer()

    peaks = [(peak.getMZ(), peak.getIntensity()) for peak in container]
    if not peaks:
        return []

    max_int = max(p[1] for p in peaks)
    if max_int > 0:
        peaks = [(mz, intensity / max_int * 100.0) for mz, intensity in peaks]

    return peaks


def score_pattern(
    observed: list[tuple[float, float]],
    theoretical: list[tuple[float, float]],
) -> dict:
    """Score observed vs theoretical isotope patterns.

    Parameters
    ----------
    observed:
        Observed (mz, intensity) pairs.
    theoretical:
        Theoretical (mz, intensity) pairs.

    Returns
    -------
    dict
        Contains cosine_score, per-peak comparisons.
    """
    n = min(len(observed), len(theoretical))
    if n == 0:
        return {"cosine_score": 0.0, "n_peaks_compared": 0, "peaks": []}

    obs_int = [observed[i][1] for i in range(n)]
    theo_int = [theoretical[i][1] for i in range(n)]

    # Normalize observed to max=100
    obs_max = max(obs_int) if obs_int else 1.0
    if obs_max > 0:
        obs_int = [v / obs_max * 100.0 for v in obs_int]

    # Cosine similarity
    dot = sum(a * b for a, b in zip(obs_int, theo_int))
    mag_a = sum(a ** 2 for a in obs_int) ** 0.5
    mag_b = sum(b ** 2 for b in theo_int) ** 0.5

    cosine = dot / (mag_a * mag_b) if (mag_a > 0 and mag_b > 0) else 0.0

    peaks = []
    for i in range(n):
        peaks.append({
            "peak_index": i,
            "obs_mz": round(observed[i][0], 6),
            "theo_mz": round(theoretical[i][0], 6),
            "obs_intensity": round(obs_int[i], 4),
            "theo_intensity": round(theo_int[i], 4),
        })

    return {
        "cosine_score": round(cosine, 6),
        "n_peaks_compared": n,
        "peaks": peaks,
    }


@click.command()
@click.option("--observed", required=True,
              help='Observed pattern as "mz:int,mz:int,..." (e.g. "180.063:100,181.067:6.5")')
@click.option("--formula", required=True, help="Molecular formula (e.g. C6H12O6)")
@click.option("--output", required=True, help="Output JSON file")
def main(observed, formula, output):
    observed_peaks = parse_observed(observed)
    theoretical = get_theoretical_pattern(formula, n_peaks=len(observed_peaks))
    result = score_pattern(observed_peaks, theoretical)
    result["formula"] = formula

    with open(output, "w") as fh:
        json.dump(result, fh, indent=2)

    print(f"Isotope fit written to {output}")
    print(f"  Cosine score: {result['cosine_score']:.6f}")
    print(f"  Peaks compared: {result['n_peaks_compared']}")


if __name__ == "__main__":
    main()
