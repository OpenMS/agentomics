"""
Blank Subtraction Tool
=======================
Subtract blank features from sample features based on intensity
fold-change thresholds and m/z+RT matching.

Features in the sample that are also present in the blank (within
tolerances) and do not exceed the fold-change threshold are removed.

Usage
-----
    python blank_subtraction_tool.py --sample sample_features.tsv --blank blank_features.tsv \
        --fold-change 3 --output cleaned.tsv
"""

import csv
import sys

import click

try:
    import pyopenms as oms  # noqa: F401
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def subtract_blanks(
    sample_features: list[dict],
    blank_features: list[dict],
    fold_change: float = 3.0,
    mz_tolerance_ppm: float = 10.0,
    rt_tolerance: float = 10.0,
) -> list[dict]:
    """Remove sample features that are present in the blank.

    Parameters
    ----------
    sample_features:
        List of dicts with keys: mz, rt, intensity.
    blank_features:
        List of dicts with keys: mz, rt, intensity.
    fold_change:
        Minimum sample/blank intensity ratio to keep a feature.
    mz_tolerance_ppm:
        m/z matching tolerance in ppm.
    rt_tolerance:
        RT matching tolerance in seconds.

    Returns
    -------
    list[dict]
        Cleaned sample features not attributable to blank.
    """
    cleaned = []

    for sf in sample_features:
        s_mz = float(sf["mz"])
        s_rt = float(sf["rt"])
        s_int = float(sf["intensity"])

        is_blank = False
        for bf in blank_features:
            b_mz = float(bf["mz"])
            b_rt = float(bf["rt"])
            b_int = float(bf["intensity"])

            mz_tol_da = s_mz * mz_tolerance_ppm / 1e6
            if abs(s_mz - b_mz) <= mz_tol_da and abs(s_rt - b_rt) <= rt_tolerance:
                if b_int > 0 and s_int / b_int < fold_change:
                    is_blank = True
                    break

        if not is_blank:
            result = dict(sf)
            result["blank_subtracted"] = "kept"
            cleaned.append(result)

    return cleaned


@click.command()
@click.option("--sample", required=True, help="Sample features TSV")
@click.option("--blank", required=True, help="Blank features TSV")
@click.option("--fold-change", type=float, default=3.0,
              help="Minimum sample/blank fold-change to keep (default: 3)")
@click.option("--mz-tolerance", type=float, default=10.0,
              help="m/z tolerance in ppm (default: 10)")
@click.option("--rt-tolerance", type=float, default=10.0,
              help="RT tolerance in seconds (default: 10)")
@click.option("--output", required=True, help="Output cleaned TSV")
def main(sample, blank, fold_change, mz_tolerance, rt_tolerance, output):
    sample_data = []
    with open(sample) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sample_data.append(row)

    blank_data = []
    with open(blank) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            blank_data.append(row)

    cleaned = subtract_blanks(
        sample_data, blank_data,
        fold_change=fold_change,
        mz_tolerance_ppm=mz_tolerance,
        rt_tolerance=rt_tolerance,
    )

    fieldnames = list(cleaned[0].keys()) if cleaned else ["mz", "rt", "intensity", "blank_subtracted"]
    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(cleaned)

    removed = len(sample_data) - len(cleaned)
    print(f"Blank subtraction: {len(sample_data)} input, {removed} removed, {len(cleaned)} kept")
    print(f"Output written to {output}")


if __name__ == "__main__":
    main()
