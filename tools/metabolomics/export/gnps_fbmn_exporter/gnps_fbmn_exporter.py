"""
GNPS FBMN Exporter
==================
Export MS2 spectra and a quantification table in GNPS Feature-Based Molecular
Networking (FBMN) format. Reads an mzML file for spectra and a feature table
(TSV) for quantification. Produces an MGF file (with SCANS=feature_id,
PEPMASS, RTINSECONDS) and a quantification CSV.

Usage
-----
    python gnps_fbmn_exporter.py --mzml data.mzML --features features.tsv \
        --output-mgf gnps.mgf --output-quant quant.csv
"""

import csv
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def load_mzml(path: str) -> oms.MSExperiment:
    """Load an mzML file into an MSExperiment object.

    Parameters
    ----------
    path:
        Path to the mzML file.

    Returns
    -------
    pyopenms.MSExperiment
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(path, exp)
    return exp


def load_features(path: str) -> list[dict]:
    """Load a feature table from TSV.

    Expected columns: feature_id, mz, rt, intensity.
    Optional additional sample columns for quantification.

    Parameters
    ----------
    path:
        Path to TSV file.

    Returns
    -------
    list of dict
    """
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            parsed = {}
            for key, val in row.items():
                try:
                    parsed[key] = float(val)
                except (ValueError, TypeError):
                    parsed[key] = val
            rows.append(parsed)
    return rows


def find_best_ms2(
    exp: oms.MSExperiment,
    mz: float,
    rt: float,
    mz_tol: float = 0.01,
    rt_tol: float = 30.0,
) -> list[tuple[float, float]] | None:
    """Find the best matching MS2 spectrum for a feature.

    Parameters
    ----------
    exp:
        MSExperiment containing spectra.
    mz:
        Precursor m/z of the feature.
    rt:
        Retention time (seconds) of the feature.
    mz_tol:
        m/z tolerance in Da for precursor matching.
    rt_tol:
        RT tolerance in seconds.

    Returns
    -------
    list of (mz, intensity) or None
        Peak list of the best matching MS2 spectrum, or None if no match.
    """
    best_spectrum = None
    best_intensity = -1.0

    for spec in exp:
        if spec.getMSLevel() != 2:
            continue
        spec_rt = spec.getRT()
        if abs(spec_rt - rt) > rt_tol:
            continue
        precursors = spec.getPrecursors()
        if not precursors:
            continue
        prec_mz = precursors[0].getMZ()
        if abs(prec_mz - mz) > mz_tol:
            continue
        # Pick the one with highest TIC
        tic = sum(spec.get_peaks()[1]) if spec.size() > 0 else 0.0
        if tic > best_intensity:
            best_intensity = tic
            mzs, intensities = spec.get_peaks()
            best_spectrum = list(zip(mzs.tolist(), intensities.tolist()))

    return best_spectrum


def write_mgf(
    features: list[dict],
    spectra: dict[str, list[tuple[float, float]]],
    path: str,
) -> int:
    """Write an MGF file in GNPS FBMN format.

    Parameters
    ----------
    features:
        Feature table rows.
    spectra:
        Mapping from feature_id to peak list [(mz, intensity), ...].
    path:
        Output MGF path.

    Returns
    -------
    int
        Number of spectra written.
    """
    count = 0
    with open(path, "w") as fh:
        for feat in features:
            fid = str(feat.get("feature_id", ""))
            if fid not in spectra:
                continue
            peaks = spectra[fid]
            if not peaks:
                continue
            mz = float(feat["mz"])
            rt = float(feat["rt"])
            fh.write("BEGIN IONS\n")
            fh.write(f"SCANS={fid}\n")
            fh.write(f"PEPMASS={mz:.6f}\n")
            fh.write(f"RTINSECONDS={rt:.2f}\n")
            fh.write("CHARGE=1+\n")
            for peak_mz, peak_int in peaks:
                fh.write(f"{peak_mz:.6f}\t{peak_int:.4f}\n")
            fh.write("END IONS\n\n")
            count += 1
    return count


def write_quant_table(features: list[dict], path: str) -> None:
    """Write a quantification table in CSV format for GNPS.

    Parameters
    ----------
    features:
        Feature table rows.
    path:
        Output CSV path.
    """
    if not features:
        with open(path, "w") as fh:
            fh.write("# No features\n")
        return
    # Ensure 'row ID' is present (GNPS expects it)
    output_rows = []
    for feat in features:
        row = dict(feat)
        if "row ID" not in row:
            row["row ID"] = row.get("feature_id", "")
        if "row m/z" not in row:
            row["row m/z"] = row.get("mz", "")
        if "row retention time" not in row:
            row["row retention time"] = row.get("rt", "")
        output_rows.append(row)

    out_fields = list(output_rows[0].keys())
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(output_rows)


def export_fbmn(
    mzml_path: str,
    features: list[dict],
    mgf_path: str,
    quant_path: str,
    mz_tol: float = 0.01,
    rt_tol: float = 30.0,
) -> tuple[int, int]:
    """Full FBMN export pipeline.

    Parameters
    ----------
    mzml_path:
        Path to mzML file.
    features:
        Feature table rows.
    mgf_path:
        Output MGF path.
    quant_path:
        Output quantification CSV path.
    mz_tol:
        m/z tolerance for MS2 matching.
    rt_tol:
        RT tolerance (seconds) for MS2 matching.

    Returns
    -------
    tuple of (spectra_written, features_count)
    """
    exp = load_mzml(mzml_path)
    spectra = {}
    for feat in features:
        fid = str(feat.get("feature_id", ""))
        mz = float(feat["mz"])
        rt = float(feat["rt"])
        ms2 = find_best_ms2(exp, mz, rt, mz_tol=mz_tol, rt_tol=rt_tol)
        if ms2 is not None:
            spectra[fid] = ms2

    n_spectra = write_mgf(features, spectra, mgf_path)
    write_quant_table(features, quant_path)
    return n_spectra, len(features)


@click.command()
@click.option("--mzml", required=True, help="Input mzML file")
@click.option("--features", required=True, help="Feature table (TSV) with feature_id, mz, rt, intensity")
@click.option("--output-mgf", required=True, help="Output MGF file")
@click.option("--output-quant", required=True, help="Output quantification CSV")
@click.option("--mz-tol", type=float, default=0.01, help="m/z tolerance in Da (default: 0.01)")
@click.option("--rt-tol", type=float, default=30.0, help="RT tolerance in seconds (default: 30)")
def main(mzml, features, output_mgf, output_quant, mz_tol, rt_tol) -> None:
    """CLI entry point."""
    features_data = load_features(features)
    n_spectra, n_features = export_fbmn(
        mzml, features_data, output_mgf, output_quant,
        mz_tol=mz_tol, rt_tol=rt_tol,
    )
    print(f"Exported {n_spectra} MS2 spectra and {n_features} features for GNPS FBMN")


if __name__ == "__main__":
    main()
