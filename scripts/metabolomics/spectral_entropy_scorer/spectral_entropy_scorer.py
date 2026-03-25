"""
Spectral Entropy Scorer
========================
Compute spectral entropy and entropy-based similarity between mass spectra.
Implements the entropy similarity score from Li & Fiehn (2021).

Usage
-----
    python spectral_entropy_scorer.py --query query_peaks.tsv --library lib_peaks.tsv \\
        --tolerance 0.02 --output scores.tsv
"""

import argparse
import csv
import math
import sys

try:
    import pyopenms as oms  # noqa: F401
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

try:
    import numpy as np  # noqa: F401
except ImportError:
    sys.exit("numpy is required. Install it with:  pip install numpy")


def normalize_intensities(intensities: list) -> list:
    """Normalize intensities to sum to 1.

    Parameters
    ----------
    intensities:
        List of intensity values.

    Returns
    -------
    list of normalized intensities (weights).
    """
    total = sum(intensities)
    if total == 0:
        return [0.0] * len(intensities)
    return [i / total for i in intensities]


def spectral_entropy(mzs: list, intensities: list) -> float:
    """Compute the spectral entropy of a spectrum.

    Entropy = -sum(w_i * log(w_i)) where w_i = intensity_i / sum(intensities).

    Parameters
    ----------
    mzs:
        List of m/z values.
    intensities:
        List of intensity values.

    Returns
    -------
    float: Spectral entropy value.
    """
    weights = normalize_intensities(intensities)
    entropy = 0.0
    for w in weights:
        if w > 0:
            entropy -= w * math.log(w)
    return entropy


def match_peaks(
    mzs_a: list,
    ints_a: list,
    mzs_b: list,
    ints_b: list,
    tolerance: float = 0.02,
) -> tuple:
    """Match peaks between two spectra within a mass tolerance.

    Parameters
    ----------
    mzs_a, ints_a:
        m/z and intensity arrays for spectrum A.
    mzs_b, ints_b:
        m/z and intensity arrays for spectrum B.
    tolerance:
        m/z tolerance in Daltons.

    Returns
    -------
    tuple of (merged_mzs, merged_ints_a, merged_ints_b)
        where unmatched peaks have 0 intensity in the other spectrum.
    """
    used_b = set()
    matched_a = []
    matched_b_idx = []

    # For each peak in A, find closest match in B
    for i, mz_a in enumerate(mzs_a):
        best_j = -1
        best_diff = tolerance + 1
        for j, mz_b in enumerate(mzs_b):
            if j in used_b:
                continue
            diff = abs(mz_a - mz_b)
            if diff <= tolerance and diff < best_diff:
                best_j = j
                best_diff = diff

        matched_a.append(i)
        if best_j >= 0:
            matched_b_idx.append(best_j)
            used_b.add(best_j)
        else:
            matched_b_idx.append(-1)

    # Build merged arrays
    merged_mzs = []
    merged_a = []
    merged_b = []

    for idx_a, idx_b in zip(matched_a, matched_b_idx):
        merged_mzs.append(mzs_a[idx_a])
        merged_a.append(ints_a[idx_a])
        merged_b.append(ints_b[idx_b] if idx_b >= 0 else 0.0)

    # Add unmatched B peaks
    for j in range(len(mzs_b)):
        if j not in used_b:
            merged_mzs.append(mzs_b[j])
            merged_a.append(0.0)
            merged_b.append(ints_b[j])

    return merged_mzs, merged_a, merged_b


def entropy_similarity(
    mzs_a: list,
    ints_a: list,
    mzs_b: list,
    ints_b: list,
    tolerance: float = 0.02,
) -> float:
    """Compute entropy similarity between two spectra (Li & Fiehn 2021).

    Entropy similarity = 1 - (2 * H_merged - H_a - H_b) / log(4)

    where H_merged is the entropy of the merged (averaged) spectrum,
    and H_a, H_b are entropies of the individual spectra.

    Parameters
    ----------
    mzs_a, ints_a:
        Query spectrum peaks.
    mzs_b, ints_b:
        Library spectrum peaks.
    tolerance:
        m/z tolerance in Daltons.

    Returns
    -------
    float: Entropy similarity score in [0, 1].
    """
    if not ints_a or not ints_b:
        return 0.0

    merged_mzs, merged_a, merged_b = match_peaks(mzs_a, ints_a, mzs_b, ints_b, tolerance)

    # Normalize each spectrum
    w_a = normalize_intensities(merged_a)
    w_b = normalize_intensities(merged_b)

    # Compute individual entropies
    h_a = 0.0
    for w in w_a:
        if w > 0:
            h_a -= w * math.log(w)

    h_b = 0.0
    for w in w_b:
        if w > 0:
            h_b -= w * math.log(w)

    # Compute merged spectrum entropy
    w_merged = [(a + b) / 2.0 for a, b in zip(w_a, w_b)]
    h_merged = 0.0
    for w in w_merged:
        if w > 0:
            h_merged -= w * math.log(w)

    # Entropy similarity
    d_ab = 2 * h_merged - h_a - h_b
    similarity = 1.0 - d_ab / math.log(4)

    # Clamp to [0, 1]
    return max(0.0, min(1.0, similarity))


def read_peaks_file(path: str) -> list:
    """Read a peaks file with columns: spectrum_id, mz, intensity.

    Parameters
    ----------
    path:
        TSV file path.

    Returns
    -------
    list of dicts with keys: spectrum_id, mzs (list), intensities (list).
    """
    spectra = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            sid = row["spectrum_id"]
            if sid not in spectra:
                spectra[sid] = {"spectrum_id": sid, "mzs": [], "intensities": []}
            spectra[sid]["mzs"].append(float(row["mz"]))
            spectra[sid]["intensities"].append(float(row["intensity"]))
    return list(spectra.values())


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Compute spectral entropy similarity between query and library spectra."
    )
    parser.add_argument("--query", required=True,
                        help="TSV with query peaks (spectrum_id, mz, intensity).")
    parser.add_argument("--library", required=True,
                        help="TSV with library peaks (spectrum_id, mz, intensity).")
    parser.add_argument("--tolerance", type=float, default=0.02,
                        help="m/z tolerance in Da (default: 0.02).")
    parser.add_argument("--output", required=True, help="Output TSV with similarity scores.")
    args = parser.parse_args()

    query_spectra = read_peaks_file(args.query)
    library_spectra = read_peaks_file(args.library)

    results = []
    for qs in query_spectra:
        q_entropy = spectral_entropy(qs["mzs"], qs["intensities"])
        for ls in library_spectra:
            score = entropy_similarity(
                qs["mzs"], qs["intensities"],
                ls["mzs"], ls["intensities"],
                tolerance=args.tolerance,
            )
            results.append({
                "query_id": qs["spectrum_id"],
                "library_id": ls["spectrum_id"],
                "query_entropy": round(q_entropy, 4),
                "entropy_similarity": round(score, 4),
            })

    fieldnames = ["query_id", "library_id", "query_entropy", "entropy_similarity"]
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    print(f"Computed {len(results)} pairwise scores, wrote to {args.output}")


if __name__ == "__main__":
    main()
