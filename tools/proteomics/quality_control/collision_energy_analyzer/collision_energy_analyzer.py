"""
Collision Energy Analyzer
=========================
Extract collision energy (CE) values from MS2 spectra in mzML files and
produce a summary report.

Usage
-----
    python collision_energy_analyzer.py --input run.mzML --output ce_analysis.tsv
"""

import csv
from typing import List

import click
import pyopenms as oms


def extract_collision_energies(exp: oms.MSExperiment) -> List[dict]:
    """Extract collision energy values from MS2 spectra.

    Parameters
    ----------
    exp:
        Loaded MSExperiment object.

    Returns
    -------
    list
        List of dicts with spectrum index, RT, precursor m/z, charge, and CE.
    """
    results = []
    for i, spec in enumerate(exp.getSpectra()):
        if spec.getMSLevel() < 2:
            continue

        rt = spec.getRT()
        precursors = spec.getPrecursors()

        for prec in precursors:
            mz = prec.getMZ()
            charge = prec.getCharge()
            ce = prec.getMetaValue("collision energy") if prec.metaValueExists("collision energy") else None

            # Also check activation energy
            if ce is None:
                activation = prec.getActivationEnergy()
                if activation > 0:
                    ce = activation

            results.append({
                "spectrum_index": i,
                "rt": round(rt, 4),
                "precursor_mz": round(mz, 6),
                "charge": charge,
                "collision_energy": round(float(ce), 2) if ce is not None else "N/A",
            })

    return results


def summarize_ce(records: List[dict]) -> dict:
    """Compute summary statistics for collision energies.

    Parameters
    ----------
    records:
        List of dicts from extract_collision_energies().

    Returns
    -------
    dict
        Summary with count, unique CEs, min, max, mean.
    """
    ce_values = [r["collision_energy"] for r in records if r["collision_energy"] != "N/A"]
    if not ce_values:
        return {"total_ms2": len(records), "with_ce": 0, "unique_ce": [], "min_ce": None, "max_ce": None}

    return {
        "total_ms2": len(records),
        "with_ce": len(ce_values),
        "unique_ce": sorted(set(ce_values)),
        "min_ce": min(ce_values),
        "max_ce": max(ce_values),
        "mean_ce": round(sum(ce_values) / len(ce_values), 2),
    }


def write_tsv(records: List[dict], output_path: str) -> None:
    """Write CE records to TSV.

    Parameters
    ----------
    records:
        List of record dicts.
    output_path:
        Output file path.
    """
    if not records:
        return
    fieldnames = list(records[0].keys())
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in records:
            writer.writerow(row)


@click.command(help="Extract collision energy values from mzML MS2 spectra.")
@click.option("--input", "input", required=True, help="Input mzML file")
@click.option("--output", default=None, help="Output TSV file path")
def main(input, output):
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input, exp)

    records = extract_collision_energies(exp)
    summary = summarize_ce(records)

    print(f"Total MS2 spectra : {summary['total_ms2']}")
    print(f"With CE values    : {summary.get('with_ce', 0)}")
    if summary.get("min_ce") is not None:
        print(f"CE range          : {summary['min_ce']} - {summary['max_ce']}")
        print(f"Mean CE           : {summary['mean_ce']}")
        print(f"Unique CE values  : {summary['unique_ce']}")

    if output:
        write_tsv(records, output)
        print(f"\nResults written to {output}")


if __name__ == "__main__":
    main()
