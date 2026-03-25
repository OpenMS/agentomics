"""
RT Prediction (Additive Model)
================================
Predict peptide retention times using additive hydrophobicity models.

Features
--------
- Krokhin model retention coefficients
- Meek model retention coefficients
- Per-residue contribution to predicted RT
- Support for modified sequences via pyopenms

Usage
-----
    python rt_prediction_additive.py --sequence PEPTIDEK --model krokhin
    python rt_prediction_additive.py --sequence PEPTIDEK --model meek --output prediction.json
"""

import csv
import json
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

# Krokhin retention coefficients (simplified version of SSRCalc)
KROKHIN_COEFFICIENTS = {
    "A": 0.62, "R": -0.60, "N": -0.60, "D": -0.46, "C": 0.29,
    "E": -0.11, "Q": -0.73, "G": 0.05, "H": -0.24, "I": 3.21,
    "L": 3.61, "K": -0.79, "M": 1.83, "F": 3.54, "P": 0.26,
    "S": -0.38, "T": -0.07, "W": 4.30, "Y": 1.31, "V": 1.87,
}

# Meek RP-HPLC retention coefficients
MEEK_COEFFICIENTS = {
    "A": 0.5, "R": -1.1, "N": -0.64, "D": -0.27, "C": -0.02,
    "E": -0.05, "Q": -0.91, "G": 0.0, "H": 0.14, "I": 2.46,
    "L": 2.46, "K": -1.54, "M": 1.31, "F": 2.65, "P": 0.38,
    "S": -0.18, "T": 0.01, "W": 3.23, "Y": 0.96, "V": 1.14,
}

MODELS = {
    "krokhin": KROKHIN_COEFFICIENTS,
    "meek": MEEK_COEFFICIENTS,
}


def predict_rt(sequence: str, model: str = "krokhin") -> dict:
    """Predict retention time using an additive model.

    Parameters
    ----------
    sequence : str
        Peptide sequence (plain or modified pyopenms notation).
    model : str
        Model name ('krokhin' or 'meek').

    Returns
    -------
    dict
        Dictionary with predicted RT, per-residue contributions, and model info.
    """
    aa_seq = oms.AASequence.fromString(sequence)
    plain = aa_seq.toUnmodifiedString()

    coefficients = MODELS.get(model, KROKHIN_COEFFICIENTS)

    contributions = []
    total = 0.0
    for i, aa in enumerate(plain):
        coeff = coefficients.get(aa, 0.0)
        total += coeff
        contributions.append({
            "position": i + 1,
            "residue": aa,
            "coefficient": coeff,
        })

    return {
        "sequence": sequence,
        "unmodified_sequence": plain,
        "model": model,
        "predicted_rt": round(total, 4),
        "length": len(plain),
        "residue_contributions": contributions,
    }


def predict_batch(sequences: list, model: str = "krokhin") -> list:
    """Predict RT for a batch of peptide sequences.

    Parameters
    ----------
    sequences : list
        List of peptide sequence strings.
    model : str
        Model name.

    Returns
    -------
    list
        List of prediction result dicts.
    """
    return [predict_rt(seq, model) for seq in sequences if seq.strip()]


@click.command(help="Predict peptide RT using additive hydrophobicity models.")
@click.option("--sequence", type=str, default=None, help="Single peptide sequence.")
@click.option("--input", "input", type=str, default=None, help="TSV file with 'sequence' column.")
@click.option(
    "--model", type=click.Choice(["krokhin", "meek"]),
    default="krokhin", help="Retention model (default: krokhin).",
)
@click.option("--output", type=str, default=None, help="Output file (.json or .tsv).")
def main(sequence, input, model, output):
    """CLI entry point."""
    if not sequence and not input:
        raise click.UsageError("Provide --sequence or --input.")

    if sequence:
        result = predict_rt(sequence, model)
        if output:
            with open(output, "w") as fh:
                json.dump(result, fh, indent=2)
        else:
            print(json.dumps(result, indent=2))
    elif input:
        sequences = []
        with open(input) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                seq = row.get("sequence", "").strip()
                if seq:
                    sequences.append(seq)
        results = predict_batch(sequences, model)
        if output:
            if output.endswith(".json"):
                with open(output, "w") as fh:
                    json.dump(results, fh, indent=2)
            else:
                with open(output, "w", newline="") as fh:
                    writer = csv.DictWriter(
                        fh, fieldnames=["sequence", "model", "predicted_rt", "length"], delimiter="\t"
                    )
                    writer.writeheader()
                    for r in results:
                        writer.writerow({k: r[k] for k in ["sequence", "model", "predicted_rt", "length"]})
            print(f"Results written to {output}")
        else:
            for r in results:
                print(f"{r['sequence']}\t{r['predicted_rt']}")


if __name__ == "__main__":
    main()
