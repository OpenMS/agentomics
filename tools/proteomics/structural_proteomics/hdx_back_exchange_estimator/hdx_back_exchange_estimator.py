"""
HDX Back-Exchange Estimator
============================
Estimate per-peptide back-exchange from fully deuterated controls.

Back-exchange is the loss of deuterium during sample handling (LC separation).
By comparing expected maximum deuteration to observed fully-deuterated mass,
the per-peptide back-exchange rate can be estimated.

Usage
-----
    python hdx_back_exchange_estimator.py --peptides peptides.tsv --fully-deuterated fd.tsv \
        --max-backexchange 40 --output report.tsv
"""

import csv
from typing import Dict, List

import click
import pyopenms as oms

DEUTERIUM_MASS_SHIFT = 1.00628  # Da difference between D and H


def count_exchangeable_amides(sequence: str) -> int:
    """Count the number of exchangeable backbone amide hydrogens.

    Exchangeable amides = len(sequence) - prolines - 2.

    Parameters
    ----------
    sequence:
        Amino acid sequence string.

    Returns
    -------
    int
        Number of exchangeable amide hydrogens.
    """
    aa = oms.AASequence.fromString(sequence)
    n = aa.size()
    proline_count = sum(1 for i in range(n) if aa.getResidue(i).getOneLetterCode() == "P")
    return max(0, n - proline_count - 2)


def get_peptide_mass(sequence: str) -> float:
    """Get the monoisotopic mass of a peptide.

    Parameters
    ----------
    sequence:
        Amino acid sequence.

    Returns
    -------
    float
        Monoisotopic mass in Da.
    """
    aa = oms.AASequence.fromString(sequence)
    return aa.getMonoWeight()


def compute_theoretical_max_deuterated_mass(sequence: str) -> float:
    """Compute the theoretical mass when all exchangeable amides are deuterated.

    Parameters
    ----------
    sequence:
        Amino acid sequence.

    Returns
    -------
    float
        Theoretical fully-deuterated mass.
    """
    base_mass = get_peptide_mass(sequence)
    exchangeable = count_exchangeable_amides(sequence)
    return base_mass + exchangeable * DEUTERIUM_MASS_SHIFT


def compute_back_exchange(
    sequence: str,
    undeuterated_mass: float,
    fully_deuterated_mass: float,
) -> Dict[str, float]:
    """Compute back-exchange for a single peptide.

    Parameters
    ----------
    sequence:
        Peptide sequence.
    undeuterated_mass:
        Observed undeuterated centroid mass.
    fully_deuterated_mass:
        Observed fully-deuterated centroid mass.

    Returns
    -------
    dict
        Dictionary with back_exchange_da, back_exchange_pct, exchangeable_amides, etc.
    """
    exchangeable = count_exchangeable_amides(sequence)
    theoretical_max = undeuterated_mass + exchangeable * DEUTERIUM_MASS_SHIFT

    observed_shift = fully_deuterated_mass - undeuterated_mass
    theoretical_shift = exchangeable * DEUTERIUM_MASS_SHIFT

    if theoretical_shift <= 0:
        back_exchange_pct = 0.0
        back_exchange_da = 0.0
    else:
        back_exchange_da = theoretical_shift - observed_shift
        back_exchange_pct = (back_exchange_da / theoretical_shift) * 100.0

    return {
        "sequence": sequence,
        "exchangeable_amides": exchangeable,
        "undeuterated_mass": round(undeuterated_mass, 6),
        "fully_deuterated_mass": round(fully_deuterated_mass, 6),
        "theoretical_max_mass": round(theoretical_max, 6),
        "observed_shift_da": round(observed_shift, 6),
        "theoretical_shift_da": round(theoretical_shift, 6),
        "back_exchange_da": round(back_exchange_da, 6),
        "back_exchange_pct": round(back_exchange_pct, 2),
    }


def flag_high_back_exchange(
    results: List[Dict[str, float]],
    max_backexchange: float = 40.0,
) -> List[Dict[str, object]]:
    """Flag peptides exceeding maximum allowed back-exchange.

    Parameters
    ----------
    results:
        List of back-exchange result dicts.
    max_backexchange:
        Maximum allowed back-exchange percentage.

    Returns
    -------
    list
        Results with added 'flag' column.
    """
    flagged = []
    for r in results:
        row = dict(r)
        row["exceeds_threshold"] = "YES" if r["back_exchange_pct"] > max_backexchange else "NO"
        flagged.append(row)
    return flagged


def read_peptides(peptides_path: str) -> Dict[str, float]:
    """Read peptides TSV with columns: sequence, centroid_mass."""
    peptides = {}
    with open(peptides_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            peptides[row["sequence"]] = float(row["centroid_mass"])
    return peptides


def read_fully_deuterated(fd_path: str) -> Dict[str, float]:
    """Read fully deuterated controls TSV with columns: sequence, centroid_mass."""
    fd = {}
    with open(fd_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            fd[row["sequence"]] = float(row["centroid_mass"])
    return fd


def write_output(output_path: str, results: List[Dict[str, object]]) -> None:
    """Write back-exchange report to TSV."""
    if not results:
        return
    fieldnames = list(results[0].keys())
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


@click.command(help="Estimate per-peptide back-exchange from fully deuterated controls.")
@click.option("--peptides", required=True, help="Undeuterated peptides TSV (sequence, centroid_mass)")
@click.option("--fully-deuterated", required=True, help="Fully deuterated TSV (sequence, centroid_mass)")
@click.option(
    "--max-backexchange", type=float, default=40.0,
    help="Maximum allowed back-exchange percentage (default: 40)",
)
@click.option("--output", required=True, help="Output report TSV file")
def main(peptides, fully_deuterated, max_backexchange, output):
    peptides_data = read_peptides(peptides)
    fd = read_fully_deuterated(fully_deuterated)

    results = []
    for seq, undeut_mass in peptides_data.items():
        if seq in fd:
            result = compute_back_exchange(seq, undeut_mass, fd[seq])
            results.append(result)

    flagged = flag_high_back_exchange(results, max_backexchange)
    write_output(output, flagged)

    n_flagged = sum(1 for r in flagged if r["exceeds_threshold"] == "YES")
    print(f"Processed {len(flagged)} peptides")
    print(f"Peptides exceeding {max_backexchange}% back-exchange: {n_flagged}")
    print(f"Output written to {output}")


if __name__ == "__main__":
    main()
