"""
Spectrum Similarity Scorer
==========================
Compute cosine similarity between two MS2 spectra from MGF files.
Uses a custom MGF reader and pyopenms SpectrumAlignment for peak matching.

Features:
- Custom MGF file parser (no MascotGenericFile dependency)
- Cosine similarity scoring with configurable tolerance
- TSV output with query_id, library_id, score, matched_peaks

Usage
-----
    python spectrum_similarity_scorer.py --query query.mgf --library reference.mgf --tolerance 0.02
    python spectrum_similarity_scorer.py --query query.mgf --library ref.mgf --tolerance 0.02 --output scores.tsv
"""

import csv
import math
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def parse_mgf(filepath: str) -> list[dict]:
    """Parse an MGF file into a list of spectrum dictionaries.

    Parameters
    ----------
    filepath : str
        Path to the MGF file.

    Returns
    -------
    list[dict]
        List of dicts with keys: title, pepmass, charge, mz_array, intensity_array.
    """
    spectra = []
    current = None

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line == "BEGIN IONS":
                current = {"title": "", "pepmass": 0.0, "charge": 0, "mz_array": [], "intensity_array": []}
            elif line == "END IONS":
                if current is not None:
                    spectra.append(current)
                    current = None
            elif current is not None:
                if line.startswith("TITLE="):
                    current["title"] = line[6:]
                elif line.startswith("PEPMASS="):
                    parts = line[8:].split()
                    current["pepmass"] = float(parts[0])
                elif line.startswith("CHARGE="):
                    charge_str = line[7:].replace("+", "").replace("-", "")
                    current["charge"] = int(charge_str) if charge_str else 0
                elif line and line[0].isdigit():
                    parts = line.split()
                    if len(parts) >= 2:
                        current["mz_array"].append(float(parts[0]))
                        current["intensity_array"].append(float(parts[1]))

    return spectra


def mgf_to_msspectrum(spec_dict: dict) -> oms.MSSpectrum:
    """Convert a parsed MGF spectrum dict to an MSSpectrum object.

    Parameters
    ----------
    spec_dict : dict
        Parsed spectrum dictionary from parse_mgf.

    Returns
    -------
    oms.MSSpectrum
        pyopenms MSSpectrum object.
    """
    spectrum = oms.MSSpectrum()
    spectrum.set_peaks((spec_dict["mz_array"], spec_dict["intensity_array"]))
    spectrum.sortByPosition()
    return spectrum


def cosine_similarity(query: oms.MSSpectrum, library: oms.MSSpectrum, tolerance: float = 0.02) -> dict:
    """Compute cosine similarity between two spectra.

    Parameters
    ----------
    query : oms.MSSpectrum
        Query spectrum.
    library : oms.MSSpectrum
        Library/reference spectrum.
    tolerance : float
        Mass tolerance in Da for peak matching (default 0.02).

    Returns
    -------
    dict
        Dictionary with keys: score, matched_peaks.
    """
    aligner = oms.SpectrumAlignment()
    param = aligner.getParameters()
    param.setValue("tolerance", tolerance)
    param.setValue("is_relative_tolerance", "false")
    aligner.setParameters(param)

    alignment = []
    aligner.getSpectrumAlignment(alignment, query, library)

    if not alignment:
        return {"score": 0.0, "matched_peaks": 0}

    q_mz, q_int = query.get_peaks()
    l_mz, l_int = library.get_peaks()

    dot_product = 0.0
    q_norm = 0.0
    l_norm = 0.0

    matched_q_indices = set()
    matched_l_indices = set()

    for qi, li in alignment:
        q_i_val = math.sqrt(q_int[qi])
        l_i_val = math.sqrt(l_int[li])
        dot_product += q_i_val * l_i_val
        matched_q_indices.add(qi)
        matched_l_indices.add(li)

    for i in range(len(q_int)):
        q_norm += q_int[i]
    for i in range(len(l_int)):
        l_norm += l_int[i]

    q_norm = math.sqrt(q_norm)
    l_norm = math.sqrt(l_norm)

    if q_norm == 0 or l_norm == 0:
        return {"score": 0.0, "matched_peaks": len(alignment)}

    score = dot_product / (q_norm * l_norm)
    return {"score": round(score, 6), "matched_peaks": len(alignment)}


def score_spectra(
    query_path: str, library_path: str, tolerance: float = 0.02
) -> list[dict]:
    """Score all query spectra against all library spectra.

    Parameters
    ----------
    query_path : str
        Path to query MGF file.
    library_path : str
        Path to library MGF file.
    tolerance : float
        Mass tolerance in Da for peak matching.

    Returns
    -------
    list[dict]
        List of dicts with keys: query_id, library_id, score, matched_peaks.
    """
    query_spectra = parse_mgf(query_path)
    library_spectra = parse_mgf(library_path)

    results = []
    for qi, q_dict in enumerate(query_spectra):
        q_spec = mgf_to_msspectrum(q_dict)
        q_id = q_dict["title"] if q_dict["title"] else f"query_{qi}"
        for li, l_dict in enumerate(library_spectra):
            l_spec = mgf_to_msspectrum(l_dict)
            l_id = l_dict["title"] if l_dict["title"] else f"library_{li}"
            sim = cosine_similarity(q_spec, l_spec, tolerance)
            results.append({
                "query_id": q_id,
                "library_id": l_id,
                "score": sim["score"],
                "matched_peaks": sim["matched_peaks"],
            })

    return results


def write_tsv(results: list[dict], output_path: str) -> None:
    """Write scoring results to TSV file.

    Parameters
    ----------
    results : list[dict]
        List of scoring result dictionaries.
    output_path : str
        Path to output TSV file.
    """
    fieldnames = ["query_id", "library_id", "score", "matched_peaks"]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


@click.command(help="Compute cosine similarity between MS2 spectra from MGF files.")
@click.option("--query", required=True, help="Path to query MGF file")
@click.option("--library", required=True, help="Path to library/reference MGF file")
@click.option("--tolerance", type=float, default=0.02, help="Mass tolerance in Da (default: 0.02)")
@click.option("--output", default=None, help="Output TSV file path (default: print to stdout)")
def main(query, library, tolerance, output):
    results = score_spectra(query, library, tolerance)

    if output:
        write_tsv(results, output)
        print(f"Wrote {len(results)} scores to {output}")
    else:
        print("query_id\tlibrary_id\tscore\tmatched_peaks")
        for r in results:
            print(f"{r['query_id']}\t{r['library_id']}\t{r['score']}\t{r['matched_peaks']}")


if __name__ == "__main__":
    main()
