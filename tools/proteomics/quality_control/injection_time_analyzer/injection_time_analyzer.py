"""
Injection Time Analyzer
========================
Extract ion injection time values from mzML spectrum metadata.

Injection times are stored as spectrum-level metadata (CV parameter
MS:1000927 "ion injection time"). This tool extracts and summarizes them.

Usage
-----
    python injection_time_analyzer.py --input run.mzML --output injection_times.tsv
"""

import csv
import math
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

ION_INJECTION_TIME_ACCESSION = "MS:1000927"


def extract_injection_times(exp: oms.MSExperiment) -> list[dict]:
    """Extract injection times from all spectra in an experiment.

    Parameters
    ----------
    exp:
        Loaded ``pyopenms.MSExperiment``.

    Returns
    -------
    list[dict]
        Each dict has: scan_index, rt, ms_level, injection_time_ms.
    """
    results = []
    for i, spec in enumerate(exp.getSpectra()):
        injection_time = None

        # Try to get from float data arrays (common storage)
        if spec.metaValueExists(ION_INJECTION_TIME_ACCESSION):
            injection_time = float(spec.getMetaValue(ION_INJECTION_TIME_ACCESSION))
        elif spec.metaValueExists("ion injection time"):
            injection_time = float(spec.getMetaValue("ion injection time"))

        # Also check instrument settings
        if injection_time is None:
            acq = spec.getAcquisitionInfo()
            if acq:
                for a in acq:
                    if a.metaValueExists(ION_INJECTION_TIME_ACCESSION):
                        injection_time = float(a.getMetaValue(ION_INJECTION_TIME_ACCESSION))
                        break

        results.append({
            "scan_index": i,
            "rt": round(spec.getRT(), 4),
            "ms_level": spec.getMSLevel(),
            "injection_time_ms": round(injection_time, 4) if injection_time is not None else None,
        })

    return results


def summarize_injection_times(records: list[dict]) -> dict:
    """Summarize injection times by MS level.

    Parameters
    ----------
    records:
        Output of ``extract_injection_times``.

    Returns
    -------
    dict
        Per-level statistics.
    """
    by_level: dict[int, list[float]] = {}
    for r in records:
        if r["injection_time_ms"] is not None:
            level = r["ms_level"]
            by_level.setdefault(level, []).append(r["injection_time_ms"])

    summary = {}
    for level, times in sorted(by_level.items()):
        mean = sum(times) / len(times)
        std = math.sqrt(sum((t - mean) ** 2 for t in times) / len(times)) if times else 0.0
        summary[f"MS{level}"] = {
            "count": len(times),
            "mean_ms": round(mean, 4),
            "std_ms": round(std, 4),
            "min_ms": round(min(times), 4),
            "max_ms": round(max(times), 4),
        }

    return summary


@click.command(help="Extract injection time values from mzML metadata.")
@click.option("--input", "input", required=True, help="Path to mzML file")
@click.option("--output", required=True, help="Output TSV file")
def main(input, output):
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input, exp)

    records = extract_injection_times(exp)

    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["scan_index", "rt", "ms_level", "injection_time_ms"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(records)

    n_with_time = sum(1 for r in records if r["injection_time_ms"] is not None)
    print(f"Wrote {len(records)} scans to {output} ({n_with_time} with injection times)")

    summary = summarize_injection_times(records)
    for level, stats in summary.items():
        print(f"  {level}: mean={stats['mean_ms']:.1f} ms, std={stats['std_ms']:.1f} ms")


if __name__ == "__main__":
    main()
