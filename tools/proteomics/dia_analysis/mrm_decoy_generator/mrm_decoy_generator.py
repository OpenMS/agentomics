"""
MRM Decoy Generator
===================
Generate decoy transitions from target transition lists using DecoyGenerator.

Features
--------
- Shuffle or reverse peptide sequences to create decoys
- Preserve precursor/product relationships
- Prefix decoy entries with configurable tag

Usage
-----
    python mrm_decoy_generator.py --input transitions.tsv --output decoy_transitions.tsv --method shuffle
    python mrm_decoy_generator.py --input transitions.tsv --output decoy_transitions.tsv --method reverse
"""

import csv
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

PROTON = 1.007276


def _generate_decoy_sequence(sequence: str, method: str, seed: int = 42) -> str:
    """Generate a decoy peptide sequence using DecoyGenerator.

    Parameters
    ----------
    sequence : str
        Target peptide sequence.
    method : str
        Method: 'shuffle' or 'reverse'.
    seed : int
        Random seed for shuffle.

    Returns
    -------
    str
        Decoy peptide sequence.
    """
    dg = oms.DecoyGenerator()
    dg.setSeed(seed)
    aa = oms.AASequence.fromString(sequence)

    if method == "reverse":
        decoy_aa = dg.reversePeptides(aa, b"Trypsin")
        return decoy_aa.toUnmodifiedString()
    else:
        results = dg.shuffle(aa, b"Trypsin", 1)
        if results:
            return results[0].toUnmodifiedString()
        return sequence[::-1]  # fallback


def _compute_product_mzs(sequence: str, charge: int, n_fragments: int = 6) -> list:
    """Compute theoretical product m/z values for a peptide.

    Parameters
    ----------
    sequence : str
        Peptide sequence.
    charge : int
        Precursor charge state.
    n_fragments : int
        Maximum number of y-ion fragments to return.

    Returns
    -------
    list
        List of product m/z values (y-ions, charge 1).
    """
    aa_seq = oms.AASequence.fromString(sequence)

    tsg = oms.TheoreticalSpectrumGenerator()
    params = tsg.getParameters()
    params.setValue("add_y_ions", "true")
    params.setValue("add_b_ions", "false")
    params.setValue("add_metainfo", "true")
    tsg.setParameters(params)

    spec = oms.MSSpectrum()
    tsg.getSpectrum(spec, aa_seq, 1, charge)

    mzs = []
    for i in range(min(spec.size(), n_fragments)):
        mzs.append(spec[i].getMZ())

    return mzs


def generate_mrm_decoys(
    input_path: str,
    output_path: str,
    method: str = "shuffle",
    decoy_tag: str = "DECOY_",
) -> int:
    """Generate decoy transitions from a target transition list.

    Reads target transitions from TSV, generates decoy peptide sequences,
    computes theoretical product m/z values for decoys, and writes the
    decoy transitions to a new TSV.

    Parameters
    ----------
    input_path : str
        Path to input transitions TSV.
    output_path : str
        Path to output decoy transitions TSV.
    method : str
        Decoy method: 'shuffle' or 'reverse'.
    decoy_tag : str
        Prefix tag for decoy identifiers.

    Returns
    -------
    int
        Number of decoy transitions generated.
    """
    # Read target transitions
    targets = []
    with open(input_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames
        for row in reader:
            targets.append(row)

    # Group targets by transition_group_id (or PeptideSequence)
    groups = {}
    for t in targets:
        gid = t.get("transition_group_id", t.get("PeptideSequence", ""))
        if gid not in groups:
            groups[gid] = []
        groups[gid].append(t)

    decoys = []
    seed = 42

    for group_id, group_trans in groups.items():
        sequence = group_trans[0].get("PeptideSequence", "")
        if not sequence:
            continue

        decoy_seq = _generate_decoy_sequence(sequence, method, seed)
        seed += 1

        # Compute precursor m/z for decoy
        try:
            aa_decoy = oms.AASequence.fromString(decoy_seq)
            charge = int(group_trans[0].get("PrecursorCharge", 2))
            decoy_prec_mz = (aa_decoy.getMonoWeight() + charge * PROTON) / charge
        except Exception:
            decoy_prec_mz = float(group_trans[0].get("PrecursorMz", 0.0))

        # Compute product m/z values for decoy
        decoy_product_mzs = _compute_product_mzs(decoy_seq, charge, len(group_trans))

        for i, target in enumerate(group_trans):
            decoy_row = dict(target)
            decoy_row["PeptideSequence"] = decoy_seq
            decoy_row["PrecursorMz"] = str(round(decoy_prec_mz, 6))

            if i < len(decoy_product_mzs):
                decoy_row["ProductMz"] = str(round(decoy_product_mzs[i], 6))

            if "ProteinName" in decoy_row:
                decoy_row["ProteinName"] = decoy_tag + decoy_row["ProteinName"]
            if "transition_name" in decoy_row:
                decoy_row["transition_name"] = decoy_tag + decoy_row["transition_name"]
            if "transition_group_id" in decoy_row:
                decoy_row["transition_group_id"] = (
                    decoy_tag + decoy_row["transition_group_id"]
                )

            decoys.append(decoy_row)

    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(decoys)

    return len(decoys)


@click.command(help="Generate MRM decoy transitions from a target transition list.")
@click.option("--input", "input_path", required=True, help="Path to target transitions TSV.")
@click.option("--output", required=True, help="Path to output decoy transitions TSV.")
@click.option(
    "--method", type=click.Choice(["shuffle", "reverse"]), default="shuffle",
    help="Decoy generation method (default: shuffle).",
)
@click.option("--decoy-tag", default="DECOY_", help="Decoy tag prefix (default: DECOY_).")
def main(input_path, output, method, decoy_tag):
    """CLI entry point."""
    n = generate_mrm_decoys(input_path, output, method, decoy_tag)
    print(f"Generated {n} decoy transitions -> {output}")


if __name__ == "__main__":
    main()
