"""
SCP Reporter QC
================
Quality control for single-cell proteomics (SCP) data using isobaric
reporter ions.  Computes the sample-to-carrier ratio per spectrum, which
is a key QC metric for carrier-based SCP experiments.

The carrier channel typically has much higher intensity than single-cell
channels.  Ratios that are too high indicate insufficient carrier signal;
ratios that are too low suggest excessive carrier relative to single cells.

Usage
-----
    python scp_reporter_qc.py --input reporter_ions.tsv \
        --carrier-channel 131C --output qc.tsv
"""

import argparse
import csv
import math
import sys
from typing import Dict, List

try:
    import pyopenms as oms  # noqa: F401
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def compute_sample_to_carrier_ratios(
    spectra: List[Dict[str, float]], carrier_channel: str
) -> List[Dict[str, object]]:
    """Compute sample-to-carrier ratio for each spectrum.

    Parameters
    ----------
    spectra:
        List of dicts mapping channel name to intensity.  Must include
        a ``spectrum_id`` key (string).
    carrier_channel:
        Name of the carrier channel (e.g. ``"131C"``).

    Returns
    -------
    list of dict
        One entry per spectrum with ``spectrum_id``, ``carrier_intensity``,
        ``mean_sample_intensity``, ``sample_to_carrier_ratio``,
        ``num_nonzero_samples``.
    """
    results: List[Dict[str, object]] = []
    for spectrum in spectra:
        spec_id = spectrum.get("spectrum_id", "unknown")
        carrier_int = spectrum.get(carrier_channel, 0.0)

        sample_intensities = []
        for ch, val in spectrum.items():
            if ch in ("spectrum_id", carrier_channel):
                continue
            if isinstance(val, (int, float)) and val > 0:
                sample_intensities.append(val)

        mean_sample = sum(sample_intensities) / len(sample_intensities) if sample_intensities else 0.0
        ratio = mean_sample / carrier_int if carrier_int > 0 else float("nan")

        results.append({
            "spectrum_id": spec_id,
            "carrier_intensity": carrier_int,
            "mean_sample_intensity": mean_sample,
            "sample_to_carrier_ratio": ratio,
            "num_nonzero_samples": len(sample_intensities),
        })
    return results


def qc_summary(ratios: List[Dict[str, object]]) -> Dict[str, object]:
    """Compute summary statistics over sample-to-carrier ratios.

    Returns
    -------
    dict
        ``n_spectra``, ``median_ratio``, ``mean_ratio``, ``std_ratio``,
        ``below_0_01_count`` (spectra with ratio < 0.01, possibly problematic).
    """
    valid_ratios = [
        r["sample_to_carrier_ratio"] for r in ratios
        if isinstance(r["sample_to_carrier_ratio"], float)
        and not math.isnan(r["sample_to_carrier_ratio"])
    ]
    if not valid_ratios:
        return {
            "n_spectra": len(ratios),
            "median_ratio": float("nan"),
            "mean_ratio": float("nan"),
            "std_ratio": float("nan"),
            "below_0_01_count": 0,
        }

    valid_ratios_sorted = sorted(valid_ratios)
    n = len(valid_ratios_sorted)
    median = valid_ratios_sorted[n // 2] if n % 2 == 1 else (
        (valid_ratios_sorted[n // 2 - 1] + valid_ratios_sorted[n // 2]) / 2.0
    )
    mean = sum(valid_ratios) / n
    variance = sum((r - mean) ** 2 for r in valid_ratios) / n if n > 1 else 0.0
    std = math.sqrt(variance)

    below_threshold = sum(1 for r in valid_ratios if r < 0.01)

    return {
        "n_spectra": len(ratios),
        "median_ratio": median,
        "mean_ratio": mean,
        "std_ratio": std,
        "below_0_01_count": below_threshold,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Single-cell proteomics QC: sample-to-carrier ratio per spectrum."
    )
    parser.add_argument(
        "--input", required=True,
        help="Input TSV with spectrum_id and reporter ion intensities per channel",
    )
    parser.add_argument(
        "--carrier-channel", required=True,
        help="Name of the carrier channel (e.g. 131C)",
    )
    parser.add_argument("--output", required=True, help="Output QC TSV")
    args = parser.parse_args()

    spectra: List[Dict[str, float]] = []
    with open(args.input, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            spec: Dict[str, float] = {}
            for key, val in row.items():
                if key == "spectrum_id":
                    spec["spectrum_id"] = val
                else:
                    try:
                        spec[key] = float(val)
                    except (ValueError, TypeError):
                        spec[key] = 0.0
            spectra.append(spec)

    if not spectra:
        sys.exit("No spectra found in input.")

    if args.carrier_channel not in (spectra[0] if spectra else {}):
        print(f"Warning: carrier channel '{args.carrier_channel}' not found in input columns.")

    ratios = compute_sample_to_carrier_ratios(spectra, args.carrier_channel)
    summary = qc_summary(ratios)

    with open(args.output, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "spectrum_id", "carrier_intensity", "mean_sample_intensity",
            "sample_to_carrier_ratio", "num_nonzero_samples",
        ])
        for r in ratios:
            ratio_str = f"{r['sample_to_carrier_ratio']:.6f}" if not math.isnan(
                r["sample_to_carrier_ratio"]
            ) else "NA"
            writer.writerow([
                r["spectrum_id"],
                f"{r['carrier_intensity']:.2f}",
                f"{r['mean_sample_intensity']:.2f}",
                ratio_str,
                r["num_nonzero_samples"],
            ])
        writer.writerow([])
        writer.writerow(["metric", "value"])
        for key, val in summary.items():
            if isinstance(val, float) and not math.isnan(val):
                writer.writerow([key, f"{val:.6f}"])
            else:
                writer.writerow([key, val])

    median_str = f"{summary['median_ratio']:.4f}" if not math.isnan(summary["median_ratio"]) else "NA"
    print(f"Processed {summary['n_spectra']} spectra, median ratio: {median_str}")


if __name__ == "__main__":
    main()
