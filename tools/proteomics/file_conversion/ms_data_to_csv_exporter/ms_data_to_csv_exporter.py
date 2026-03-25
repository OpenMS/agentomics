"""
MS Data to CSV Exporter
=======================
Export mzML or featureXML data to flat CSV/TSV files.

Usage
-----
    python ms_data_to_csv_exporter.py --input data.mzML --type peaks --output peaks.tsv
    python ms_data_to_csv_exporter.py --input features.featureXML --type features --output features.tsv
"""

import csv
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def export_mzml_peaks(input_path: str, output_path: str, ms_level: int = 0) -> dict:
    """Export peak data from an mzML file to TSV.

    Parameters
    ----------
    input_path : str
        Path to input mzML file.
    output_path : str
        Path to output TSV file.
    ms_level : int
        MS level to filter (0 = all levels).

    Returns
    -------
    dict
        Export statistics.
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    total_peaks = 0
    spectra_count = 0

    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["spectrum_index", "ms_level", "rt", "mz", "intensity"])

        for i, spectrum in enumerate(exp):
            level = spectrum.getMSLevel()
            if ms_level > 0 and level != ms_level:
                continue
            rt = spectrum.getRT()
            mz_array, intensity_array = spectrum.get_peaks()
            spectra_count += 1
            for mz, intensity in zip(mz_array, intensity_array):
                writer.writerow([i, level, f"{rt:.4f}", f"{mz:.6f}", f"{intensity:.4f}"])
                total_peaks += 1

    return {"spectra_exported": spectra_count, "total_peaks": total_peaks}


def export_mzml_spectra_summary(input_path: str, output_path: str) -> dict:
    """Export spectrum-level summary from an mzML file to TSV."""
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["index", "native_id", "ms_level", "rt", "num_peaks", "base_peak_mz", "tic"])

        for i, spectrum in enumerate(exp):
            mz_array, intensity_array = spectrum.get_peaks()
            base_peak_mz = 0.0
            tic = 0.0
            if len(intensity_array) > 0:
                max_idx = intensity_array.argmax()
                base_peak_mz = mz_array[max_idx]
                tic = float(intensity_array.sum())

            writer.writerow([
                i, spectrum.getNativeID(), spectrum.getMSLevel(),
                f"{spectrum.getRT():.4f}", len(mz_array),
                f"{base_peak_mz:.6f}", f"{tic:.4f}"
            ])

    return {"spectra_exported": exp.size()}


def export_featurexml(input_path: str, output_path: str) -> dict:
    """Export features from a featureXML file to TSV."""
    feature_map = oms.FeatureMap()
    oms.FeatureXMLFile().load(input_path, feature_map)

    with open(output_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["feature_id", "rt", "mz", "intensity", "charge", "quality"])

        for feature in feature_map:
            writer.writerow([
                feature.getUniqueId(),
                f"{feature.getRT():.4f}",
                f"{feature.getMZ():.6f}",
                f"{feature.getIntensity():.4f}",
                feature.getCharge(),
                f"{feature.getOverallQuality():.4f}",
            ])

    return {"features_exported": feature_map.size()}


@click.command(help="Export mzML or featureXML data to flat TSV.")
@click.option("--input", "input", required=True, help="Input file (mzML or featureXML)")
@click.option(
    "--type", "type", required=True, type=click.Choice(["peaks", "spectra", "features"]),
    help="Export type: peaks (mzML), spectra (mzML summary), features (featureXML)",
)
@click.option("--ms-level", type=int, default=0, help="MS level filter for peaks (0=all)")
@click.option("--output", required=True, help="Output TSV file")
def main(input, type, ms_level, output) -> None:
    if type == "peaks":
        stats = export_mzml_peaks(input, output, ms_level)
        print(f"Exported {stats['total_peaks']} peaks from {stats['spectra_exported']} spectra")
    elif type == "spectra":
        stats = export_mzml_spectra_summary(input, output)
        print(f"Exported {stats['spectra_exported']} spectra summaries")
    elif type == "features":
        stats = export_featurexml(input, output)
        print(f"Exported {stats['features_exported']} features")


if __name__ == "__main__":
    main()
