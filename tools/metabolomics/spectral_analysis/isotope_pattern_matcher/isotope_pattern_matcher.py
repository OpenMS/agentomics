"""
Isotope Pattern Generator & Matcher
=====================================
Generate the theoretical isotope distribution for a molecular formula
and optionally compare it against an observed spectrum to compute a
cosine similarity score using pyopenms.

Usage
-----
    # Generate isotope pattern for glucose
    python isotope_pattern_matcher.py --formula C6H12O6

    # Compare against observed peaks (m/z intensity pairs on stdin or --peaks)
    python isotope_pattern_matcher.py --formula C6H12O6 \\
        --peaks 181.0709,100.0 --peaks 182.0742,6.7 --peaks 183.0775,0.4
"""

import math
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit(
        "pyopenms is required. Install it with:  pip install pyopenms"
    )


def get_isotope_distribution(
    formula: str,
    max_isotopes: int = 5,
) -> list:
    """Compute the theoretical isotope distribution for a molecular formula.

    Parameters
    ----------
    formula:
        Empirical formula string, e.g. ``"C6H12O6"``.
    max_isotopes:
        Maximum number of isotope peaks to return (default: 5).

    Returns
    -------
    list of (mz, relative_abundance)
        Sorted by m/z; relative abundances sum to 100.
    """
    ef = oms.EmpiricalFormula(formula)
    isotope_dist = ef.getIsotopeDistribution(
        oms.CoarseIsotopePatternGenerator(max_isotopes)
    )
    peaks = [(p.getMZ(), p.getIntensity()) for p in isotope_dist.getContainer()]
    if not peaks:
        return []
    max_ab = max(ab for _, ab in peaks)
    return [(mz, ab / max_ab * 100) for mz, ab in peaks]


def cosine_similarity(
    theoretical: list,
    observed: list,
    mz_tolerance: float = 0.02,
) -> float:
    """Compute cosine similarity between theoretical and observed peak lists.

    Parameters
    ----------
    theoretical:
        List of ``(mz, intensity)`` tuples for the theoretical pattern.
    observed:
        List of ``(mz, intensity)`` tuples for the observed spectrum.
    mz_tolerance:
        Maximum m/z difference to match a pair of peaks (default: 0.02 Da).

    Returns
    -------
    float
        Cosine similarity in [0, 1].
    """
    theo_vec = []
    obs_vec = []
    for tmz, tint in theoretical:
        matched = 0.0
        for omz, oint in observed:
            if abs(omz - tmz) <= mz_tolerance:
                matched = oint
                break
        theo_vec.append(tint)
        obs_vec.append(matched)

    dot = sum(t * o for t, o in zip(theo_vec, obs_vec))
    norm_t = math.sqrt(sum(t ** 2 for t in theo_vec))
    norm_o = math.sqrt(sum(o ** 2 for o in obs_vec))
    if norm_t == 0 or norm_o == 0:
        return 0.0
    return dot / (norm_t * norm_o)


def parse_peaks(peak_strings: list) -> list:
    """Parse ``"mz,intensity"`` strings into ``(float, float)`` tuples."""
    peaks = []
    for s in peak_strings:
        parts = s.split(",")
        if len(parts) != 2:
            raise ValueError(
                f"Expected 'mz,intensity' format, got: {s!r}"
            )
        peaks.append((float(parts[0]), float(parts[1])))
    return peaks


@click.command()
@click.option("--formula", required=True, help="Molecular formula (e.g. C6H12O6)")
@click.option("--max-isotopes", type=int, default=5,
              help="Maximum isotope peaks to compute (default: 5)")
@click.option("--peaks", multiple=True,
              help="Observed peaks as 'mz,intensity' pairs for similarity scoring")
@click.option("--tolerance", type=float, default=0.02,
              help="m/z tolerance in Da for peak matching (default: 0.02)")
def main(formula, max_isotopes, peaks, tolerance):
    distribution = get_isotope_distribution(formula, max_isotopes)
    if not distribution:
        print("Could not compute isotope distribution for the given formula.")
        return

    print(f"Isotope distribution for {formula}:")
    print(f"\n{'Peak':>5}  {'m/z':>12}  {'Relative Abundance (%)':>22}")
    print("-" * 44)
    for i, (mz, rel_ab) in enumerate(distribution):
        bar = "#" * int(rel_ab / 5)
        print(f"  M+{i}  {mz:>12.4f}  {rel_ab:>6.2f} %  {bar}")

    if peaks:
        observed = parse_peaks(list(peaks))
        sim = cosine_similarity(distribution, observed, tolerance)
        print(f"\nCosine similarity vs. observed peaks: {sim:.4f}")
        if sim >= 0.9:
            print("  ✓ Excellent match (≥ 0.90)")
        elif sim >= 0.7:
            print("  ~ Good match (≥ 0.70)")
        else:
            print("  ✗ Poor match (< 0.70)")


if __name__ == "__main__":
    main()
