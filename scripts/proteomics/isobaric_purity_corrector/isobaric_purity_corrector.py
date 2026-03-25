"""
Isobaric Purity Corrector
==========================
Correct TMT or iTRAQ reporter ion quantification for isotopic impurity
using a manufacturer-provided purity correction matrix.  The tool reads
a quantification TSV and a purity matrix CSV, solves the linear system
to produce corrected intensities.

The correction is: corrected = inv(purity_matrix) @ observed

Usage
-----
    python isobaric_purity_corrector.py --input quant.tsv \
        --label TMT16plex --purity-matrix purity.csv --output corrected.tsv
"""

import argparse
import csv
import sys
from typing import Dict, List, Tuple

try:
    import pyopenms as oms  # noqa: F401
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

import numpy as np

# Default channel names for common labeling schemes
LABEL_CHANNELS: Dict[str, List[str]] = {
    "TMT6plex": ["126", "127", "128", "129", "130", "131"],
    "TMT10plex": ["126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131"],
    "TMT11plex": [
        "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C",
    ],
    "TMT16plex": [
        "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N",
        "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N",
    ],
    "TMT18plex": [
        "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N",
        "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N", "134C", "135N",
    ],
    "iTRAQ4plex": ["114", "115", "116", "117"],
    "iTRAQ8plex": ["113", "114", "115", "116", "117", "118", "119", "121"],
}


def load_purity_matrix(purity_path: str) -> np.ndarray:
    """Load purity correction matrix from CSV.

    The CSV should have N rows and N columns (no headers), where N is the
    number of channels.  Each row represents the contribution of a single
    channel to all observed channels.

    Parameters
    ----------
    purity_path:
        Path to purity matrix CSV.

    Returns
    -------
    numpy.ndarray
        Square purity matrix of shape (N, N).
    """
    rows: List[List[float]] = []
    with open(purity_path, newline="") as fh:
        reader = csv.reader(fh)
        for row in reader:
            values = [float(v.strip()) for v in row if v.strip()]
            if values:
                rows.append(values)
    return np.array(rows, dtype=float)


def correct_intensities(
    observed: np.ndarray, purity_matrix: np.ndarray
) -> np.ndarray:
    """Apply purity correction to observed intensities.

    Parameters
    ----------
    observed:
        Array of shape (n_spectra, n_channels) with observed intensities.
    purity_matrix:
        Square purity matrix of shape (n_channels, n_channels).

    Returns
    -------
    numpy.ndarray
        Corrected intensities, same shape as observed.  Negative values
        are clipped to zero.
    """
    # Solve: purity_matrix @ corrected = observed (for each spectrum)
    inv_matrix = np.linalg.inv(purity_matrix)
    corrected = observed @ inv_matrix.T
    corrected[corrected < 0] = 0.0
    return corrected


def get_channels(label: str) -> List[str]:
    """Return channel names for a given labeling scheme.

    Parameters
    ----------
    label:
        Labeling scheme name (e.g. ``"TMT16plex"``).

    Returns
    -------
    list of str
        Channel names.
    """
    if label not in LABEL_CHANNELS:
        raise ValueError(
            f"Unknown label '{label}'. Supported: {', '.join(sorted(LABEL_CHANNELS.keys()))}"
        )
    return LABEL_CHANNELS[label]


def process_quant_file(
    input_path: str,
    purity_matrix: np.ndarray,
    channels: List[str],
) -> Tuple[List[str], List[Dict[str, str]], np.ndarray]:
    """Read quantification file and apply purity correction.

    Parameters
    ----------
    input_path:
        Path to input TSV.
    purity_matrix:
        Purity correction matrix.
    channels:
        Channel names to correct.

    Returns
    -------
    tuple
        (non_channel_fields, metadata_rows, corrected_array)
    """
    metadata_rows: List[Dict[str, str]] = []
    observed_list: List[List[float]] = []

    with open(input_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        all_fields = reader.fieldnames or []
        non_channel_fields = [f for f in all_fields if f not in channels]

        for row in reader:
            meta = {f: row.get(f, "") for f in non_channel_fields}
            metadata_rows.append(meta)
            intensities = []
            for ch in channels:
                val = row.get(ch, "0").strip()
                try:
                    intensities.append(float(val))
                except (ValueError, TypeError):
                    intensities.append(0.0)
            observed_list.append(intensities)

    observed = np.array(observed_list, dtype=float)
    corrected = correct_intensities(observed, purity_matrix)
    return non_channel_fields, metadata_rows, corrected


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Correct TMT/iTRAQ reporter ion quantification for isotopic impurity."
    )
    parser.add_argument("--input", required=True, help="Input quantification TSV")
    parser.add_argument(
        "--label", required=True,
        help="Labeling scheme (e.g. TMT6plex, TMT10plex, TMT16plex, iTRAQ4plex)",
    )
    parser.add_argument("--purity-matrix", required=True, help="Purity correction matrix CSV")
    parser.add_argument("--output", required=True, help="Output corrected TSV")
    args = parser.parse_args()

    channels = get_channels(args.label)
    purity_matrix = load_purity_matrix(args.purity_matrix)

    expected_size = len(channels)
    if purity_matrix.shape != (expected_size, expected_size):
        sys.exit(
            f"Purity matrix shape {purity_matrix.shape} does not match "
            f"{expected_size} channels for {args.label}"
        )

    non_ch_fields, meta_rows, corrected = process_quant_file(args.input, purity_matrix, channels)

    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(non_ch_fields + channels)
        for i, meta in enumerate(meta_rows):
            row = [meta.get(f, "") for f in non_ch_fields]
            row += [f"{corrected[i, j]:.4f}" for j in range(len(channels))]
            writer.writerow(row)

    print(f"Corrected {len(meta_rows)} spectra across {len(channels)} channels -> {args.output}")


if __name__ == "__main__":
    main()
