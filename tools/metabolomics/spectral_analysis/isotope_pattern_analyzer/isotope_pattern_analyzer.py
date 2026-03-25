"""
Isotope Pattern Analyzer
========================
Generate theoretical isotope distributions for molecular formulas, score observed
isotope patterns against theoretical ones using cosine similarity, and detect
halogenation (Cl/Br) from M+2 peak enhancement.

This tool consolidates the functionality of isotope_pattern_matcher,
isotope_pattern_scorer, and isotope_pattern_fit_scorer into a single,
more capable utility.

Features
--------
- Theoretical isotope pattern generation via pyopenms CoarseIsotopePatternGenerator
- Cosine similarity scoring between observed and theoretical patterns
- Da or ppm m/z tolerance modes
- Halogen (Cl/Br) detection from M+2 peak enhancement
- JSON output with per-peak detail
- Bar-chart preview in the terminal

Usage
-----
    # Generate isotope pattern for glucose
    python isotope_pattern_analyzer.py --formula C6H12O6

    # Score observed peaks (mz:intensity format) against formula
    python isotope_pattern_analyzer.py --formula C6H12O6 \\
        --observed "180.063:100,181.067:6.5,182.070:0.5" \\
        --output result.json

    # Score using comma-separated mz,intensity pairs (legacy format)
    python isotope_pattern_analyzer.py --formula C6H12O6 \\
        --peaks 180.063,100.0 --peaks 181.067,6.5 \\
        --output result.json

    # Use ppm tolerance
    python isotope_pattern_analyzer.py --formula C6H12O6 \\
        --observed "180.063:100,181.067:6.5" --tolerance 10 --tolerance-unit ppm
"""

import json
import math
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

try:
    import numpy as np

    _HAS_NUMPY = True
except ImportError:
    _HAS_NUMPY = False


def get_theoretical_pattern(formula: str, max_isotopes: int = 6) -> list[tuple[float, float]]:
    """Compute the theoretical isotope distribution for a molecular formula.

    Parameters
    ----------
    formula:
        Empirical formula string, e.g. ``"C6H12O6"``.
    max_isotopes:
        Maximum number of isotope peaks to return (default: 6).

    Returns
    -------
    list of (mz, relative_abundance)
        Sorted by m/z; relative abundances normalized so the maximum is 100.
    """
    ef = oms.EmpiricalFormula(formula)
    iso_dist = ef.getIsotopeDistribution(oms.CoarseIsotopePatternGenerator(max_isotopes))
    container = iso_dist.getContainer()
    peaks = [(p.getMZ(), p.getIntensity()) for p in container]
    if not peaks:
        return []
    max_int = max(p[1] for p in peaks)
    if max_int == 0:
        return peaks
    return [(mz, intensity / max_int * 100.0) for mz, intensity in peaks]


def _mz_tolerance_da(mz: float, tolerance: float, unit: str) -> float:
    """Convert tolerance to Da, supporting 'da' and 'ppm' units.

    Parameters
    ----------
    mz:
        Reference m/z value (only used when unit is 'ppm').
    tolerance:
        Tolerance value.
    unit:
        Either ``'da'`` or ``'ppm'``.

    Returns
    -------
    float
        Tolerance in Daltons.
    """
    if unit == "ppm":
        return mz * tolerance * 1e-6
    return tolerance


def cosine_similarity(
    observed: list[tuple[float, float]],
    theoretical: list[tuple[float, float]],
    tolerance: float = 0.05,
    tolerance_unit: str = "da",
) -> float:
    """Compute cosine similarity between observed and theoretical isotope patterns.

    Each theoretical peak is matched to the closest observed peak within the
    specified tolerance.  Unmatched theoretical peaks contribute zero to the
    dot product.

    Parameters
    ----------
    observed:
        List of (mz, intensity) for observed peaks.
    theoretical:
        List of (mz, intensity) for theoretical peaks.
    tolerance:
        m/z tolerance value.
    tolerance_unit:
        ``'da'`` for Daltons (default) or ``'ppm'`` for parts per million.

    Returns
    -------
    float
        Cosine similarity score in [0, 1].
    """
    if not observed or not theoretical:
        return 0.0

    if _HAS_NUMPY:
        obs_mz = np.array([p[0] for p in observed])
        obs_int = np.array([p[1] for p in observed], dtype=float)
        theo_int = np.array([p[1] for p in theoretical], dtype=float)

        matched_obs = np.zeros(len(theoretical))
        for i, (tmz, _) in enumerate(theoretical):
            tol_da = _mz_tolerance_da(tmz, tolerance, tolerance_unit)
            diffs = np.abs(obs_mz - tmz)
            mask = diffs <= tol_da
            if mask.any():
                matched_obs[i] = obs_int[np.argmin(diffs)]

        dot = float(np.dot(theo_int, matched_obs))
        norm_t = float(np.linalg.norm(theo_int))
        norm_o = float(np.linalg.norm(matched_obs))
    else:
        matched_obs = []
        for tmz, _ in theoretical:
            tol_da = _mz_tolerance_da(tmz, tolerance, tolerance_unit)
            best = 0.0
            for omz, oint in observed:
                if abs(omz - tmz) <= tol_da:
                    best = oint
                    break
            matched_obs.append(best)

        dot = sum(t * o for (_, t), o in zip(theoretical, matched_obs))
        norm_t = math.sqrt(sum(t ** 2 for _, t in theoretical))
        norm_o = math.sqrt(sum(o ** 2 for o in matched_obs))

    if norm_t == 0 or norm_o == 0:
        return 0.0
    return dot / (norm_t * norm_o)


def detect_halogenation(
    observed: list[tuple[float, float]],
    theoretical: list[tuple[float, float]],
) -> dict:
    """Detect Cl/Br halogenation from an enhanced M+2 isotope peak.

    Chlorine has a natural Cl-37 / Cl-35 ratio of ~32.5%, giving an M+2
    peak ~32.5% of M.  Bromine has an almost 1:1 Br-79/Br-81 ratio, giving
    an M+2 peak ~97% of M.  If the observed M+2 significantly exceeds the
    theoretical value, a halogen is flagged.

    Parameters
    ----------
    observed:
        Observed isotope pattern as (mz, intensity) list.
    theoretical:
        Theoretical isotope pattern (without explicit halogens) as (mz,
        intensity) list.

    Returns
    -------
    dict
        Keys: m2_ratio_observed, m2_ratio_theoretical, m2_excess,
        halogen_flag, possible_halogen.
    """
    result: dict = {
        "m2_ratio_observed": None,
        "m2_ratio_theoretical": None,
        "m2_excess": None,
        "halogen_flag": False,
        "possible_halogen": "none",
    }
    if len(observed) < 3 or len(theoretical) < 3:
        return result

    obs_m0 = observed[0][1]
    obs_m2 = observed[2][1]
    theo_m0 = max(theoretical[0][1], 1e-10)
    theo_m2 = theoretical[2][1]

    if obs_m0 == 0:
        return result

    obs_ratio = obs_m2 / obs_m0 * 100.0
    theo_ratio = theo_m2 / theo_m0 * 100.0
    excess = obs_ratio - theo_ratio

    result["m2_ratio_observed"] = round(obs_ratio, 2)
    result["m2_ratio_theoretical"] = round(theo_ratio, 2)
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


def parse_observed(observed_str: str) -> list[tuple[float, float]]:
    """Parse an observed peak string into (mz, intensity) pairs.

    Accepts both ``"mz:intensity"`` (colon-separated) and
    ``"mz,intensity"`` (comma-separated, one pair per call) formats.
    Multiple peaks are separated by commas in colon format.

    Parameters
    ----------
    observed_str:
        String of the form ``"180.063:100,181.067:6.5"`` (colon-separated
        pairs) or a single ``"180.063,100.0"`` (comma-separated pair).

    Returns
    -------
    list of (mz, intensity)

    Raises
    ------
    ValueError
        If a token cannot be parsed.
    """
    peaks = []
    # If colons are present, assume "mz:intensity,mz:intensity,..." format
    if ":" in observed_str:
        for token in observed_str.split(","):
            token = token.strip()
            if not token:
                continue
            parts = token.split(":")
            if len(parts) != 2:
                raise ValueError(f"Expected 'mz:intensity', got: {token!r}")
            peaks.append((float(parts[0].strip()), float(parts[1].strip())))
    else:
        # Single "mz,intensity" pair (legacy --peaks flag style)
        parts = observed_str.split(",")
        if len(parts) != 2:
            raise ValueError(f"Expected 'mz,intensity', got: {observed_str!r}")
        peaks.append((float(parts[0].strip()), float(parts[1].strip())))
    return peaks


def score_pattern(
    observed: list[tuple[float, float]],
    formula: str,
    max_isotopes: int = 6,
    tolerance: float = 0.05,
    tolerance_unit: str = "da",
) -> dict:
    """Score an observed isotope pattern against a theoretical one.

    Parameters
    ----------
    observed:
        List of (mz, intensity) observed peaks.
    formula:
        Molecular formula string.
    max_isotopes:
        Maximum number of theoretical isotope peaks to generate.
    tolerance:
        m/z tolerance value.
    tolerance_unit:
        ``'da'`` or ``'ppm'``.

    Returns
    -------
    dict
        Contains formula, cosine_similarity, per-peak detail, and
        halogen detection results.
    """
    theoretical = get_theoretical_pattern(formula, max_isotopes)
    cos_sim = cosine_similarity(observed, theoretical, tolerance, tolerance_unit)
    halogen = detect_halogenation(observed, theoretical)

    n = min(len(observed), len(theoretical))
    peak_detail = []
    for i in range(n):
        peak_detail.append({
            "peak_index": i,
            "obs_mz": round(observed[i][0], 6),
            "theo_mz": round(theoretical[i][0], 6),
            "obs_intensity": round(observed[i][1], 4),
            "theo_intensity": round(theoretical[i][1], 4),
        })

    return {
        "formula": formula,
        "cosine_similarity": round(cos_sim, 6),
        "n_peaks_compared": n,
        "tolerance": tolerance,
        "tolerance_unit": tolerance_unit,
        "peaks": peak_detail,
        "theoretical_pattern": [{"mz": round(mz, 6), "intensity": round(ab, 4)} for mz, ab in theoretical],
        "halogen_detection": halogen,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


@click.command(
    help="Generate theoretical isotope patterns and optionally score them against observed peaks."
)
@click.option("--formula", required=True, help="Molecular formula (e.g. C6H12O6)")
@click.option("--max-isotopes", type=int, default=6, help="Max isotope peaks to generate (default: 6)")
@click.option(
    "--observed", default=None,
    help="Observed peaks as 'mz:int,mz:int,...' (colon-separated pairs).",
)
@click.option(
    "--peaks", multiple=True,
    help="Observed peaks as 'mz,intensity' pairs (repeatable, legacy format).",
)
@click.option("--tolerance", type=float, default=0.05, help="m/z tolerance (default: 0.05)")
@click.option(
    "--tolerance-unit", type=click.Choice(["da", "ppm"]), default="da",
    help="Tolerance unit: 'da' (Daltons, default) or 'ppm'.",
)
@click.option("--output", default=None, help="Output JSON file (optional)")
def main(formula, max_isotopes, observed, peaks, tolerance, tolerance_unit, output):
    """CLI entry point."""
    theoretical = get_theoretical_pattern(formula, max_isotopes)
    if not theoretical:
        print("Could not compute isotope distribution for the given formula.")
        return

    # Print theoretical pattern
    print(f"Theoretical isotope pattern for {formula}:")
    print(f"\n{'Peak':>5}  {'m/z':>12}  {'Rel. Abundance (%)':>18}")
    print("-" * 42)
    for i, (mz, rel_ab) in enumerate(theoretical):
        bar = "#" * int(rel_ab / 5)
        print(f"  M+{i}  {mz:>12.4f}  {rel_ab:>6.2f} %  {bar}")

    # Collect observed peaks
    obs_peaks: list[tuple[float, float]] = []
    if observed:
        obs_peaks = parse_observed(observed)
    elif peaks:
        for p in peaks:
            obs_peaks.extend(parse_observed(p))

    if obs_peaks:
        result = score_pattern(obs_peaks, formula, max_isotopes, tolerance, tolerance_unit)
        cos_sim = result["cosine_similarity"]
        print(f"\nCosine similarity ({tolerance} {tolerance_unit}): {cos_sim:.4f}")
        if cos_sim >= 0.9:
            print("  ✓ Excellent match (≥ 0.90)")
        elif cos_sim >= 0.7:
            print("  ~ Good match (≥ 0.70)")
        else:
            print("  ✗ Poor match (< 0.70)")

        hal = result["halogen_detection"]
        if hal["halogen_flag"]:
            print(f"\nHalogen detected: {hal['possible_halogen']}")
            print(f"  M+2 observed: {hal['m2_ratio_observed']}% | theoretical: {hal['m2_ratio_theoretical']}%")

        if output:
            with open(output, "w") as fh:
                json.dump(result, fh, indent=2)
            print(f"\nDetailed results written to {output}")
    elif output:
        # No observed peaks but output requested — write theoretical pattern only
        data = {
            "formula": formula,
            "theoretical_pattern": [{"mz": round(mz, 6), "intensity": round(ab, 4)} for mz, ab in theoretical],
        }
        with open(output, "w") as fh:
            json.dump(data, fh, indent=2)
        print(f"\nTheoretical pattern written to {output}")


if __name__ == "__main__":
    main()
