"""
Identification QC Reporter
===========================
Report identification-level QC metrics from a peptide TSV file.

Metrics include peptide/PSM counts, score distribution, precursor mass
error statistics, and modification frequencies.

The input TSV must contain columns: sequence, charge, score, precursor_mz,
and optionally modifications.

Usage
-----
    python identification_qc_reporter.py --input results.tsv --output id_qc.json
"""

import csv
import json
import math
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

PROTON = 1.007276


def compute_id_qc(rows: list[dict]) -> dict:
    """Compute identification QC metrics from parsed peptide rows.

    Parameters
    ----------
    rows:
        List of dicts with keys: sequence, charge, score, precursor_mz,
        and optionally modifications.

    Returns
    -------
    dict
        QC metrics dictionary.
    """
    if not rows:
        return {"psm_count": 0, "unique_peptides": 0}

    sequences = [r["sequence"] for r in rows]
    scores = [float(r["score"]) for r in rows]
    unique_peptides = set(sequences)

    score_mean = sum(scores) / len(scores)
    score_std = math.sqrt(sum((s - score_mean) ** 2 for s in scores) / len(scores)) if scores else 0.0
    score_min = min(scores)
    score_max = max(scores)

    mass_errors = []
    for r in rows:
        try:
            seq = oms.AASequence.fromString(r["sequence"])
            theo_mass = seq.getMonoWeight()
            charge = int(r["charge"])
            obs_mz = float(r["precursor_mz"])
            theo_mz = (theo_mass + charge * PROTON) / charge
            error_ppm = (obs_mz - theo_mz) / theo_mz * 1e6
            mass_errors.append(error_ppm)
        except Exception:
            continue

    me_mean = sum(mass_errors) / len(mass_errors) if mass_errors else 0.0
    me_std = (
        math.sqrt(sum((e - me_mean) ** 2 for e in mass_errors) / len(mass_errors))
        if mass_errors
        else 0.0
    )

    mod_counts: dict[str, int] = {}
    for r in rows:
        mods = r.get("modifications", "").strip()
        if mods:
            for mod in mods.split(";"):
                mod = mod.strip()
                if mod:
                    mod_counts[mod] = mod_counts.get(mod, 0) + 1

    return {
        "psm_count": len(rows),
        "unique_peptides": len(unique_peptides),
        "score_mean": round(score_mean, 4),
        "score_std": round(score_std, 4),
        "score_min": round(score_min, 4),
        "score_max": round(score_max, 4),
        "mass_error_mean_ppm": round(me_mean, 4),
        "mass_error_std_ppm": round(me_std, 4),
        "modification_counts": mod_counts,
    }


def load_peptide_tsv(path: str) -> list[dict]:
    """Load a peptide TSV file.

    Parameters
    ----------
    path:
        Path to TSV file with columns: sequence, charge, score, precursor_mz.

    Returns
    -------
    list[dict]
    """
    rows = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


@click.command(help="Report identification-level QC from peptide TSV.")
@click.option("--input", "input", required=True, help="Peptide TSV file")
@click.option("--output", required=True, help="Output JSON report")
def main(input, output):
    rows = load_peptide_tsv(input)
    metrics = compute_id_qc(rows)

    with open(output, "w") as fh:
        json.dump(metrics, fh, indent=2)

    print(f"ID QC report written to {output}")
    print(f"  PSMs            : {metrics['psm_count']}")
    print(f"  Unique peptides : {metrics['unique_peptides']}")


if __name__ == "__main__":
    main()
