"""
Charge State Predictor
========================
Predict charge state distribution for peptides based on basic residues and ionization.

Features
--------
- Predict likely charge states based on number of basic residues
- Estimate charge state probabilities using Henderson-Hasselbalch
- Account for N-terminal amine and side-chain basicities
- Support for different pH conditions (ESI, MALDI)

Usage
-----
    python charge_state_predictor.py --sequence PEPTIDEK --ph 2.0 --output charges.json
"""

import argparse
import json
import sys

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

PROTON = 1.007276

# pKa values for protonatable groups in ESI conditions
PKA_VALUES = {
    "nterm": 7.7,  # alpha-amino group
    "K": 10.5,
    "R": 12.5,
    "H": 6.0,
}


def count_basic_sites(sequence: str) -> dict:
    """Count protonatable basic sites in a peptide.

    Parameters
    ----------
    sequence : str
        Amino acid sequence.

    Returns
    -------
    dict
        Counts of each type of basic site.
    """
    aa_seq = oms.AASequence.fromString(sequence)
    plain = aa_seq.toUnmodifiedString()

    sites = {"nterm": 1, "K": 0, "R": 0, "H": 0}
    for aa in plain:
        if aa in sites:
            sites[aa] += 1

    sites["total"] = sum(sites.values())
    return sites


def predict_charge_states(sequence: str, ph: float = 2.0, max_charge: int = 0) -> dict:
    """Predict charge state distribution for a peptide.

    Parameters
    ----------
    sequence : str
        Peptide sequence.
    ph : float
        Solution pH (default 2.0 for typical ESI conditions).
    max_charge : int
        Maximum charge to consider (0 = auto from basic sites).

    Returns
    -------
    dict
        Dictionary with charge state probabilities and m/z values.
    """
    aa_seq = oms.AASequence.fromString(sequence)
    plain = aa_seq.toUnmodifiedString()
    mono_mass = aa_seq.getMonoWeight()

    basic_sites = count_basic_sites(plain)
    total_basic = basic_sites["total"]

    if max_charge == 0:
        max_charge = min(total_basic, 8)  # Cap at 8

    # Calculate protonation probability for each site
    protonation_probs = {}
    protonation_probs["nterm"] = 1.0 / (1.0 + 10 ** (ph - PKA_VALUES["nterm"]))
    for aa in ["K", "R", "H"]:
        if basic_sites[aa] > 0:
            protonation_probs[aa] = 1.0 / (1.0 + 10 ** (ph - PKA_VALUES[aa]))

    # Expected charge = sum of protonation probabilities across all sites
    expected_charge = protonation_probs["nterm"]
    for aa in ["K", "R", "H"]:
        expected_charge += protonation_probs.get(aa, 0.0) * basic_sites[aa]

    # Generate charge state distribution (simplified model)
    # Use a distribution centered around expected charge
    charge_states = []
    raw_weights = []
    for z in range(1, max_charge + 1):
        # Simple Gaussian-like weighting around expected charge
        weight = 2.71828 ** (-0.5 * ((z - expected_charge) ** 2))
        raw_weights.append((z, weight))

    total_weight = sum(w for _, w in raw_weights)

    for z, weight in raw_weights:
        probability = weight / total_weight if total_weight > 0 else 0.0
        mz = (mono_mass + z * PROTON) / z
        charge_states.append({
            "charge": z,
            "probability": round(probability, 4),
            "mz": round(mz, 6),
        })

    # Sort by probability descending
    charge_states.sort(key=lambda x: x["probability"], reverse=True)

    return {
        "sequence": sequence,
        "unmodified_sequence": plain,
        "monoisotopic_mass": round(mono_mass, 6),
        "ph": ph,
        "basic_sites": basic_sites,
        "expected_charge": round(expected_charge, 2),
        "most_likely_charge": charge_states[0]["charge"] if charge_states else 1,
        "charge_states": charge_states,
    }


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(description="Predict peptide charge state distribution.")
    parser.add_argument("--sequence", required=True, help="Peptide sequence.")
    parser.add_argument("--ph", type=float, default=2.0, help="Solution pH (default: 2.0 for ESI).")
    parser.add_argument("--max-charge", type=int, default=0, help="Max charge state (0 = auto).")
    parser.add_argument("--output", type=str, help="Output file (.json).")
    args = parser.parse_args()

    result = predict_charge_states(args.sequence, args.ph, args.max_charge)

    if args.output:
        with open(args.output, "w") as fh:
            json.dump(result, fh, indent=2)
        print(f"Results written to {args.output}")
    else:
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
