"""
MassQL Query Tool
==================
Query mzML data using a simplified MassQL-like syntax.

Supported query patterns:
- ``MS2PROD=<mz>`` : Find MS2 spectra containing a product ion at the given m/z
- ``MS1MZ=<mz>`` : Find MS1 spectra containing a peak at the given m/z
- ``PRECMZ=<mz>`` : Find MS2 spectra with a specific precursor m/z

All matches use a configurable Da tolerance (default 0.5 Da).

Usage
-----
    python massql_query_tool.py --input data.mzML --query "MS2PROD=226.18" --output results.tsv
"""

import csv
import re
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def parse_query(query: str) -> dict:
    """Parse a MassQL-like query string.

    Parameters
    ----------
    query:
        Query string such as ``"MS2PROD=226.18"`` or ``"MS1MZ=180.06"``.

    Returns
    -------
    dict
        Parsed query with keys: query_type, target_mz.
    """
    query = query.strip()
    pattern = re.compile(r"^(MS2PROD|MS1MZ|PRECMZ)\s*=\s*([0-9.]+)$", re.IGNORECASE)
    match = pattern.match(query)
    if not match:
        raise ValueError(f"Unsupported query syntax: {query}")

    return {
        "query_type": match.group(1).upper(),
        "target_mz": float(match.group(2)),
    }


def execute_query(
    exp: oms.MSExperiment,
    query: dict,
    tolerance_da: float = 0.5,
) -> list[dict]:
    """Execute a parsed query against an MSExperiment.

    Parameters
    ----------
    exp:
        Loaded ``pyopenms.MSExperiment``.
    query:
        Parsed query from ``parse_query``.
    tolerance_da:
        Absolute tolerance in Da for matching.

    Returns
    -------
    list[dict]
        Matching results with scan info.
    """
    qtype = query["query_type"]
    target = query["target_mz"]
    lo = target - tolerance_da
    hi = target + tolerance_da

    results = []

    for i, spec in enumerate(exp.getSpectra()):
        level = spec.getMSLevel()

        if qtype == "MS2PROD" and level == 2:
            mzs, intensities = spec.get_peaks()
            for mz, intensity in zip(mzs, intensities):
                if lo <= mz <= hi:
                    prec_mz = spec.getPrecursors()[0].getMZ() if spec.getPrecursors() else 0.0
                    results.append({
                        "scan_index": i,
                        "rt": round(spec.getRT(), 4),
                        "ms_level": level,
                        "precursor_mz": round(prec_mz, 6),
                        "matched_mz": round(float(mz), 6),
                        "matched_intensity": round(float(intensity), 2),
                    })
                    break  # one match per spectrum

        elif qtype == "MS1MZ" and level == 1:
            mzs, intensities = spec.get_peaks()
            for mz, intensity in zip(mzs, intensities):
                if lo <= mz <= hi:
                    results.append({
                        "scan_index": i,
                        "rt": round(spec.getRT(), 4),
                        "ms_level": level,
                        "precursor_mz": 0.0,
                        "matched_mz": round(float(mz), 6),
                        "matched_intensity": round(float(intensity), 2),
                    })
                    break

        elif qtype == "PRECMZ" and level == 2:
            for prec in spec.getPrecursors():
                if lo <= prec.getMZ() <= hi:
                    results.append({
                        "scan_index": i,
                        "rt": round(spec.getRT(), 4),
                        "ms_level": level,
                        "precursor_mz": round(prec.getMZ(), 6),
                        "matched_mz": round(prec.getMZ(), 6),
                        "matched_intensity": 0.0,
                    })
                    break

    return results


@click.command()
@click.option("--input", "input_file", required=True, help="mzML file")
@click.option("--query", required=True,
              help='MassQL query (e.g. "MS2PROD=226.18", "MS1MZ=180.06", "PRECMZ=500.0")')
@click.option("--tolerance", type=float, default=0.5,
              help="m/z tolerance in Da (default: 0.5)")
@click.option("--output", required=True, help="Output results TSV")
def main(input_file, query, tolerance, output):
    parsed = parse_query(query)

    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_file, exp)

    results = execute_query(exp, parsed, tolerance_da=tolerance)

    with open(output, "w", newline="") as fh:
        fieldnames = ["scan_index", "rt", "ms_level", "precursor_mz", "matched_mz", "matched_intensity"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    print(f"Query '{query}': {len(results)} matches, written to {output}")


if __name__ == "__main__":
    main()
