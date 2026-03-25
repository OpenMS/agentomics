"""
Targeted Feature Extractor
============================
Extract chromatographic features for known target compounds from MS1
data in an mzML file.

For each target compound (defined by name, formula, and expected RT),
the tool extracts the ion chromatogram and integrates the peak area.

Usage
-----
    python targeted_feature_extractor.py --input sample.mzML --targets compounds.tsv \
        --ppm 5 --output quantified.tsv
"""

import csv

import click
import pyopenms as oms

PROTON = 1.007276


def extract_eic(
    exp: oms.MSExperiment,
    target_mz: float,
    ppm: float = 5.0,
) -> list[tuple[float, float]]:
    """Extract an extracted ion chromatogram (EIC) for a target m/z.

    Parameters
    ----------
    exp:
        Loaded ``pyopenms.MSExperiment``.
    target_mz:
        Target m/z value.
    ppm:
        Mass tolerance in ppm.

    Returns
    -------
    list[tuple[float, float]]
        List of (rt, intensity) tuples.
    """
    tol_da = target_mz * ppm / 1e6
    lo = target_mz - tol_da
    hi = target_mz + tol_da

    eic = []
    for spec in exp.getSpectra():
        if spec.getMSLevel() != 1:
            continue
        rt = spec.getRT()
        mzs, intensities = spec.get_peaks()
        total = 0.0
        for mz, intensity in zip(mzs, intensities):
            if lo <= mz <= hi:
                total += float(intensity)
        eic.append((rt, total))

    return eic


def integrate_peak(eic: list[tuple[float, float]], rt_min: float = 0.0, rt_max: float = 1e9) -> float:
    """Integrate EIC area using the trapezoidal rule.

    Parameters
    ----------
    eic:
        Extracted ion chromatogram as (rt, intensity) pairs.
    rt_min, rt_max:
        RT bounds for integration (seconds).

    Returns
    -------
    float
        Integrated area.
    """
    filtered = [(rt, i) for rt, i in eic if rt_min <= rt <= rt_max]
    if len(filtered) < 2:
        return sum(i for _, i in filtered)

    area = 0.0
    for idx in range(1, len(filtered)):
        dt = filtered[idx][0] - filtered[idx - 1][0]
        avg_int = (filtered[idx][1] + filtered[idx - 1][1]) / 2.0
        area += dt * avg_int
    return area


def extract_targets(
    exp: oms.MSExperiment,
    targets: list[dict],
    ppm: float = 5.0,
) -> list[dict]:
    """Extract features for all target compounds.

    Parameters
    ----------
    exp:
        Loaded ``pyopenms.MSExperiment``.
    targets:
        List of dicts with keys: name, formula (and optionally rt_min, rt_max).
    ppm:
        Mass tolerance in ppm.

    Returns
    -------
    list[dict]
        Each dict has: name, formula, target_mz, peak_area, max_intensity, n_points.
    """
    results = []
    for t in targets:
        formula = t["formula"]
        ef = oms.EmpiricalFormula(formula)
        neutral_mass = ef.getMonoWeight()
        target_mz = neutral_mass + PROTON  # [M+H]+

        rt_min = float(t.get("rt_min", 0.0))
        rt_max = float(t.get("rt_max", 1e9))

        eic = extract_eic(exp, target_mz, ppm=ppm)
        area = integrate_peak(eic, rt_min=rt_min, rt_max=rt_max)
        max_int = max((i for _, i in eic), default=0.0)

        results.append({
            "name": t.get("name", formula),
            "formula": formula,
            "target_mz": round(target_mz, 6),
            "peak_area": round(area, 2),
            "max_intensity": round(max_int, 2),
            "n_points": len(eic),
        })

    return results


@click.command()
@click.option("--input", "input_file", required=True, help="mzML file")
@click.option("--targets", required=True, help="Compounds TSV (name, formula)")
@click.option("--ppm", type=float, default=5.0, help="Mass tolerance in ppm (default: 5)")
@click.option("--output", required=True, help="Output quantified TSV")
def main(input_file, targets, ppm, output):
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_file, exp)

    target_list = []
    with open(targets) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            target_list.append(row)

    results = extract_targets(exp, target_list, ppm=ppm)

    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["name", "formula", "target_mz", "peak_area", "max_intensity", "n_points"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(results)

    print(f"Extracted {len(results)} targets, written to {output}")


if __name__ == "__main__":
    main()
