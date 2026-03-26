"""
FLASHDeconv Wrapper
===================
Deconvolve intact-protein mass spectra using the FLASHDeconv algorithm.
Takes an mzML file with multiply-charged spectra and produces a TSV of
deconvolved monoisotopic masses.

Wraps pyopenms.FLASHDeconvAlgorithm for basic intact mass deconvolution.

Usage
-----
    python flash_deconv.py --input intact.mzML --output masses.tsv
    python flash_deconv.py --input intact.mzML --output masses.tsv --min-mass 5000 --max-mass 50000
"""

import csv
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def deconvolve_intact(
    input_path: str,
    output_path: str,
    min_mass: float = 5000,
    max_mass: float = 100000,
) -> int:
    """Deconvolve intact protein spectra using FLASHDeconv.

    Parameters
    ----------
    input_path : str
        Path to input mzML file with multiply-charged spectra.
    output_path : str
        Path to output TSV file with deconvolved masses.
    min_mass : float
        Minimum mass to report (Da).
    max_mass : float
        Maximum mass to report (Da).

    Returns
    -------
    int
        Number of deconvolved masses found.
    """
    exp = oms.MSExperiment()
    print(f"Loading {input_path} ...")
    oms.MzMLFile().load(input_path, exp)
    exp.updateRanges()

    # Configure FLASHDeconv
    fd = oms.FLASHDeconvAlgorithm()
    params = fd.getDefaults()
    params.setValue("SD:min_mass", float(min_mass))
    params.setValue("SD:max_mass", float(max_mass))
    fd.setParameters(params)

    # Run deconvolution
    deconvolved_spectra = []
    deconvolved_features = []
    fd.run(exp, deconvolved_spectra, deconvolved_features)

    # Collect masses from deconvolved spectra via PeakGroups
    masses = []
    for dspec in deconvolved_spectra:
        for pg in dspec:
            mono_mass = pg.getMonoMass()
            if min_mass <= mono_mass <= max_mass:
                masses.append({
                    "mass": round(float(mono_mass), 4),
                    "intensity": round(float(pg.getIntensity()), 2),
                    "charge": int(pg.getRepAbsCharge()),
                    "qscore": round(float(pg.getQscore()), 4),
                })

    # Also collect from mass features
    for feat in deconvolved_features:
        mono_mass = float(feat.avg_mass)
        if min_mass <= mono_mass <= max_mass:
            masses.append({
                "mass": round(mono_mass, 4),
                "intensity": 0.0,
                "charge": 0,
                "qscore": 0.0,
            })

    # Write TSV output
    fieldnames = ["mass", "intensity", "charge", "qscore"]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(masses)

    print(f"Found {len(masses)} deconvolved masses -> {output_path}")
    return len(masses)


@click.command(help="Deconvolve intact protein spectra using FLASHDeconv.")
@click.option("--input", "input_file", required=True, help="Input mzML file")
@click.option("--output", required=True, help="Output TSV file")
@click.option(
    "--min-mass",
    type=float,
    default=5000,
    help="Minimum mass in Da (default: 5000)",
)
@click.option(
    "--max-mass",
    type=float,
    default=100000,
    help="Maximum mass in Da (default: 100000)",
)
def main(input_file, output, min_mass, max_mass):
    n = deconvolve_intact(input_file, output, min_mass, max_mass)
    print(f"Total masses reported: {n}")


if __name__ == "__main__":
    main()
