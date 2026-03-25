"""
Kovats Retention Index Calculator
=================================
Calculate Kovats retention indices from n-alkane standard retention times
for GC-MS data. Uses the logarithmic interpolation formula:

    RI = 100*n + 100 * (log(RT_x) - log(RT_n)) / (log(RT_{n+1}) - log(RT_n))

where n is the carbon number of the preceding alkane and RT_x is the
retention time of the analyte.

Usage
-----
    python kovats_ri_calculator.py --input features.tsv --standards alkane_rts.tsv --output ri_values.tsv
"""

import argparse
import csv
import math
import sys

try:
    import pyopenms as oms  # noqa: F401
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def load_tsv(path: str) -> list[dict]:
    """Load a TSV file into a list of dicts with numeric parsing.

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


def build_alkane_table(standards: list[dict]) -> list[tuple[int, float]]:
    """Build a sorted list of (carbon_number, rt) from alkane standards.

    Parameters
    ----------
    standards:
        List of dicts with keys: carbon_number, rt.

    Returns
    -------
    list of (int, float)
        Sorted by retention time.
    """
    table = []
    for row in standards:
        cn = int(row["carbon_number"])
        rt = float(row["rt"])
        if rt <= 0:
            continue
        table.append((cn, rt))
    table.sort(key=lambda x: x[1])
    return table


def calculate_kovats_ri(
    rt: float,
    alkane_table: list[tuple[int, float]],
) -> float | None:
    """Calculate the Kovats retention index for a given RT.

    Parameters
    ----------
    rt:
        Retention time of the analyte.
    alkane_table:
        Sorted list of (carbon_number, rt) for n-alkane standards.

    Returns
    -------
    float or None
        Kovats RI, or None if the RT falls outside the alkane range.
    """
    if rt <= 0:
        return None
    if len(alkane_table) < 2:
        return None

    # Find the bracketing alkanes
    for i in range(len(alkane_table) - 1):
        cn_n, rt_n = alkane_table[i]
        cn_n1, rt_n1 = alkane_table[i + 1]
        if rt_n <= rt <= rt_n1:
            if rt_n <= 0 or rt_n1 <= 0:
                return None
            log_rt = math.log10(rt)
            log_rt_n = math.log10(rt_n)
            log_rt_n1 = math.log10(rt_n1)
            denom = log_rt_n1 - log_rt_n
            if denom == 0:
                return None
            ri = 100 * cn_n + 100 * (log_rt - log_rt_n) / denom
            return round(ri, 2)
    return None


def calculate_ri_batch(
    features: list[dict],
    alkane_table: list[tuple[int, float]],
    rt_column: str = "rt",
) -> list[dict]:
    """Calculate Kovats RI for a batch of features.

    Parameters
    ----------
    features:
        List of feature dicts.
    alkane_table:
        Sorted list of (carbon_number, rt) from alkane standards.
    rt_column:
        Name of the RT column in features.

    Returns
    -------
    list of dict
        Input rows augmented with kovats_ri.
    """
    results = []
    for feat in features:
        result = dict(feat)
        rt = float(feat[rt_column])
        ri = calculate_kovats_ri(rt, alkane_table)
        result["kovats_ri"] = ri if ri is not None else "N/A"
        results.append(result)
    return results


def write_results(results: list[dict], path: str) -> None:
    """Write results to a TSV file.

    Parameters
    ----------
    results:
        List of result dicts.
    path:
        Output TSV path.
    """
    if not results:
        with open(path, "w") as fh:
            fh.write("# No results\n")
        return
    fieldnames = list(results[0].keys())
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Calculate Kovats retention indices from alkane standards for GC-MS."
    )
    parser.add_argument("--input", required=True, help="Feature table (TSV) with rt column")
    parser.add_argument("--standards", required=True, help="Alkane standards (TSV) with carbon_number, rt")
    parser.add_argument("--output", required=True, help="Output RI values (TSV)")
    parser.add_argument("--rt-column", default="rt", help="Name of RT column (default: rt)")
    args = parser.parse_args()

    features = load_tsv(args.input)
    standards = load_tsv(args.standards)
    alkane_table = build_alkane_table(standards)
    results = calculate_ri_batch(features, alkane_table, rt_column=args.rt_column)
    write_results(results, args.output)
    print(f"Calculated RI for {len(results)} features, written to {args.output}")


if __name__ == "__main__":
    main()
