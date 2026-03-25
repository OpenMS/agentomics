"""
Isotope Pattern Fit Scorer
==========================
Score observed vs theoretical isotope patterns using cosine similarity.
Detect halogenation (Cl/Br) from enhanced M+2 peaks.

Usage
-----
    python isotope_pattern_fit_scorer.py \
        --observed "180.063:100,181.067:6.5,182.070:0.5" \
        --formula C6H12O6 --output fit.json
"""

import argparse
import json
import math
import sys

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def get_theoretical_pattern(
    formula: str,
    max_isotopes: int = 6,
) -> list[tuple[float, float]]:
    """Compute the theoretical isotope distribution for a molecular formula.

    Parameters
    ----------
    formula:
        Empirical formula string, e.g. ``"C6H12O6"``.
    max_isotopes:
        Maximum number of isotope peaks.

    Returns
    -------
    list of (mz, relative_abundance)
        Relative abundances normalized so the max is 100.
    """
    ef = oms.EmpiricalFormula(formula)
    iso_dist = ef.getIsotopeDistribution(
        oms.CoarseIsotopePatternGenerator(max_isotopes)
    )
    container = iso_dist.getContainer()
    peaks = [(p.getMZ(), p.getIntensity()) for p in container]
    if not peaks:
        return []
    max_int = max(p[1] for p in peaks)
    if max_int == 0:
        return peaks
    return [(mz, intensity / max_int * 100.0) for mz, intensity in peaks]


def parse_observed(observed_str: str) -> list[tuple[float, float]]:
    """Parse observed peaks from a string.

    Format: ``"mz1:int1,mz2:int2,..."``

    Parameters
    ----------
    observed_str:
        Comma-separated mz:intensity pairs.

    Returns
    -------
    list of (mz, intensity)
    """
    peaks = []
    for pair in observed_str.split(","):
        pair = pair.strip()
        if not pair:
            continue
        parts = pair.split(":")
        if len(parts) != 2:
            raise ValueError(f"Invalid peak format: {pair!r}. Expected mz:intensity")
        mz = float(parts[0].strip())
        intensity = float(parts[1].strip())
        peaks.append((mz, intensity))
    return peaks


def cosine_similarity(
    observed: list[tuple[float, float]],
    theoretical: list[tuple[float, float]],
    mz_tolerance: float = 0.05,
) -> float:
    """Compute cosine similarity between observed and theoretical patterns.

    Parameters
    ----------
    observed:
        List of (mz, intensity) for observed peaks.
    theoretical:
        List of (mz, intensity) for theoretical peaks.
    mz_tolerance:
        Maximum m/z difference to consider a match.

    Returns
    -------
    float
        Cosine similarity score between 0 and 1.
    """
    if not observed or not theoretical:
        return 0.0

    dot_product = 0.0
    norm_obs = 0.0
    norm_theo = 0.0

    for obs_mz, obs_int in observed:
        norm_obs += obs_int ** 2

    for theo_mz, theo_int in theoretical:
        norm_theo += theo_int ** 2

    for obs_mz, obs_int in observed:
        for theo_mz, theo_int in theoretical:
            if abs(obs_mz - theo_mz) <= mz_tolerance:
                dot_product += obs_int * theo_int

    if norm_obs == 0 or norm_theo == 0:
        return 0.0
    return dot_product / (math.sqrt(norm_obs) * math.sqrt(norm_theo))


def detect_halogenation(
    observed: list[tuple[float, float]],
    theoretical: list[tuple[float, float]],
    mz_tolerance: float = 0.05,
) -> dict:
    """Detect Cl/Br halogenation from enhanced M+2 peak.

    Cl has a natural isotope ratio of ~32.5% for M+2 relative to M+0.
    Br has a natural isotope ratio of ~97% for M+2 relative to M+0.
    If the observed M+2 is significantly higher than theoretical, flag
    potential halogenation.

    Parameters
    ----------
    observed:
        Observed isotope pattern.
    theoretical:
        Theoretical isotope pattern (without halogens).
    mz_tolerance:
        m/z tolerance for peak matching.

    Returns
    -------
    dict
        Keys: m2_observed, m2_theoretical, m2_excess, halogen_flag, possible_halogen.
    """
    result = {
        "m2_observed": None,
        "m2_theoretical": None,
        "m2_excess": None,
        "halogen_flag": False,
        "possible_halogen": "none",
    }

    if len(observed) < 3 or len(theoretical) < 3:
        return result

    # M+0 is the first peak; M+2 is the third peak
    obs_m0_int = observed[0][1]
    obs_m2_int = observed[2][1] if len(observed) > 2 else 0.0
    theo_m2_int = theoretical[2][1] if len(theoretical) > 2 else 0.0

    if obs_m0_int == 0:
        return result

    obs_m2_ratio = obs_m2_int / obs_m0_int * 100.0
    theo_m2_ratio = theo_m2_int / max(theoretical[0][1], 1e-10) * 100.0

    result["m2_observed"] = round(obs_m2_ratio, 2)
    result["m2_theoretical"] = round(theo_m2_ratio, 2)
    excess = obs_m2_ratio - theo_m2_ratio
    result["m2_excess"] = round(excess, 2)

    if excess > 10.0:
        result["halogen_flag"] = True
        if excess > 70.0:
            result["possible_halogen"] = "Br"
        elif excess > 20.0:
            result["possible_halogen"] = "Cl"
        else:
            result["possible_halogen"] = "Cl (weak signal)"

    return result


def score_pattern(
    observed_str: str,
    formula: str,
    max_isotopes: int = 6,
    mz_tolerance: float = 0.05,
) -> dict:
    """Score an observed isotope pattern against a theoretical one.

    Parameters
    ----------
    observed_str:
        Comma-separated mz:intensity pairs.
    formula:
        Molecular formula.
    max_isotopes:
        Max isotope peaks.
    mz_tolerance:
        m/z tolerance.

    Returns
    -------
    dict
        Score results including cosine_similarity, halogen detection, and peak data.
    """
    observed = parse_observed(observed_str)
    theoretical = get_theoretical_pattern(formula, max_isotopes=max_isotopes)
    cos_sim = cosine_similarity(observed, theoretical, mz_tolerance=mz_tolerance)
    halogen = detect_halogenation(observed, theoretical, mz_tolerance=mz_tolerance)

    return {
        "formula": formula,
        "cosine_similarity": round(cos_sim, 6),
        "observed_peaks": [{"mz": mz, "intensity": i} for mz, i in observed],
        "theoretical_peaks": [{"mz": round(mz, 6), "intensity": round(i, 4)} for mz, i in theoretical],
        "halogen_detection": halogen,
    }


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Score observed vs theoretical isotope patterns and detect halogenation."
    )
    parser.add_argument(
        "--observed", required=True,
        help="Observed peaks as 'mz1:int1,mz2:int2,...'"
    )
    parser.add_argument("--formula", required=True, help="Molecular formula (e.g. C6H12O6)")
    parser.add_argument("--max-isotopes", type=int, default=6, help="Max isotope peaks (default: 6)")
    parser.add_argument("--mz-tolerance", type=float, default=0.05, help="m/z tolerance (default: 0.05)")
    parser.add_argument("--output", required=True, help="Output JSON file")
    args = parser.parse_args()

    result = score_pattern(
        args.observed, args.formula,
        max_isotopes=args.max_isotopes, mz_tolerance=args.mz_tolerance,
    )

    with open(args.output, "w") as fh:
        json.dump(result, fh, indent=2)

    print(f"Cosine similarity: {result['cosine_similarity']:.4f}")
    halogen = result["halogen_detection"]
    if halogen["halogen_flag"]:
        print(f"Halogen detected: {halogen['possible_halogen']} (M+2 excess: {halogen['m2_excess']}%)")
    print(f"Results written to {args.output}")


if __name__ == "__main__":
    main()
