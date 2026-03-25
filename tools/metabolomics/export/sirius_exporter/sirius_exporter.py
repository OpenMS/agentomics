"""
SIRIUS Exporter
===============
Export features and MS2 spectra to SIRIUS .ms format for molecular formula identification.

Usage
-----
    python sirius_exporter.py --features features.tsv --mzml data.mzML --output sirius_input.ms
"""

import csv
import sys
from typing import List

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def load_mzml(input_path: str) -> oms.MSExperiment:
    """Load an mzML file."""
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)
    return exp


def load_features_tsv(features_path: str) -> List[dict]:
    """Load features from a TSV file.

    Expected columns: mz, rt, charge (optional), name (optional).
    """
    features = []
    with open(features_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            feature = {
                "mz": float(row["mz"]),
                "rt": float(row["rt"]),
                "charge": int(row.get("charge", 0)) if row.get("charge") else 0,
                "name": row.get("name", ""),
            }
            features.append(feature)
    return features


def find_ms2_spectra(
    exp: oms.MSExperiment,
    precursor_mz: float,
    rt: float,
    mz_tolerance: float = 0.01,
    rt_tolerance: float = 30.0,
) -> List[oms.MSSpectrum]:
    """Find MS2 spectra matching a precursor m/z within tolerances."""
    matches = []
    for spectrum in exp:
        if spectrum.getMSLevel() != 2:
            continue
        if abs(spectrum.getRT() - rt) > rt_tolerance:
            continue
        precursors = spectrum.getPrecursors()
        if precursors:
            prec_mz = precursors[0].getMZ()
            if abs(prec_mz - precursor_mz) <= mz_tolerance:
                matches.append(spectrum)
    return matches


def write_sirius_ms(
    features: List[dict],
    exp: oms.MSExperiment,
    output_path: str,
    mz_tolerance: float = 0.01,
    rt_tolerance: float = 30.0,
) -> dict:
    """Write features and their MS2 spectra to SIRIUS .ms format.

    Returns statistics about the export.
    """
    exported = 0
    with_ms2 = 0

    with open(output_path, "w") as fh:
        for i, feature in enumerate(features):
            name = feature["name"] if feature["name"] else f"feature_{i}"
            fh.write(f">compound {name}\n")
            fh.write(f">parentmass {feature['mz']:.6f}\n")
            if feature["charge"] != 0:
                fh.write(f">charge {feature['charge']}\n")
            fh.write(f">rt {feature['rt']:.2f}\n")

            # Find matching MS2 spectra
            ms2_spectra = find_ms2_spectra(
                exp, feature["mz"], feature["rt"], mz_tolerance, rt_tolerance
            )

            if ms2_spectra:
                with_ms2 += 1
                for spectrum in ms2_spectra:
                    # Get collision energy if available
                    precursors = spectrum.getPrecursors()
                    ce = 0.0
                    if precursors:
                        ce = precursors[0].getActivationEnergy()

                    fh.write(f"\n>ms2 {feature['mz']:.6f}")
                    if ce > 0:
                        fh.write(f" {ce:.1f}")
                    fh.write("\n")

                    mz_array, intensity_array = spectrum.get_peaks()
                    for mz, intensity in zip(mz_array, intensity_array):
                        fh.write(f"{mz:.6f} {intensity:.4f}\n")

            fh.write("\n")
            exported += 1

    return {"features_exported": exported, "features_with_ms2": with_ms2}


def export_to_sirius(
    features_path: str,
    mzml_path: str,
    output_path: str,
    mz_tolerance: float = 0.01,
    rt_tolerance: float = 30.0,
) -> dict:
    """Main export function: load features and mzML, write SIRIUS .ms file."""
    features = load_features_tsv(features_path)
    exp = load_mzml(mzml_path)
    return write_sirius_ms(features, exp, output_path, mz_tolerance, rt_tolerance)


@click.command()
@click.option("--features", required=True, help="Input features TSV (columns: mz, rt, charge, name)")
@click.option("--mzml", required=True, help="Input mzML file")
@click.option("--output", required=True, help="Output SIRIUS .ms file")
@click.option("--mz-tolerance", type=float, default=0.01, help="m/z tolerance in Da (default: 0.01)")
@click.option("--rt-tolerance", type=float, default=30.0, help="RT tolerance in seconds (default: 30)")
def main(features, mzml, output, mz_tolerance, rt_tolerance) -> None:
    stats = export_to_sirius(features, mzml, output, mz_tolerance, rt_tolerance)
    print(f"Exported {stats['features_exported']} features ({stats['features_with_ms2']} with MS2) "
          f"to {output}")


if __name__ == "__main__":
    main()
