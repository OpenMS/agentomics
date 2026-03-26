"""
Decoy Database Generator
========================
Generate a target-decoy FASTA database by reversing or shuffling protein
sequences. The output contains all original (target) entries followed by
the generated decoy entries with a ``DECOY_`` identifier prefix.

Usage
-----
    python decoy_database_generator.py --input target.fasta --output target_decoy.fasta --method reverse
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit(
        "pyopenms is required. Install it with:  pip install pyopenms"
    )


def generate_decoys(input_path: str, output_path: str, method: str = "reverse") -> int:
    """Generate decoy protein sequences and write a combined target-decoy FASTA.

    Parameters
    ----------
    input_path:
        Path to the input target FASTA file.
    output_path:
        Path to the output target-decoy FASTA file.
    method:
        Decoy generation method: ``"reverse"`` or ``"shuffle"``.

    Returns
    -------
    int
        Number of decoy proteins generated.
    """
    import os

    proteins = []
    # FASTAFile.load raises on empty files, so check first
    if os.path.getsize(input_path) == 0:
        oms.FASTAFile().store(output_path, [])
        return 0

    oms.FASTAFile().load(input_path, proteins)

    if not proteins:
        oms.FASTAFile().store(output_path, [])
        return 0

    dg = oms.DecoyGenerator()
    decoys = []

    for entry in proteins:
        decoy = oms.FASTAEntry()
        decoy.identifier = "DECOY_" + entry.identifier
        decoy.description = entry.description
        seq = oms.AASequence.fromString(entry.sequence)

        if method == "reverse":
            rev = dg.reverseProtein(seq)
            decoy.sequence = rev.toString()
        elif method == "shuffle":
            shuffled_seqs = dg.shuffle(seq, "Trypsin", 1)
            # shuffle returns a list of AASequence; take the first one
            if shuffled_seqs:
                decoy.sequence = shuffled_seqs[0].toString()
            else:
                # Fallback to reverse if shuffle fails
                rev = dg.reverseProtein(seq)
                decoy.sequence = rev.toString()
        else:
            raise ValueError(f"Unknown method: {method}. Use 'reverse' or 'shuffle'.")

        decoys.append(decoy)

    combined = list(proteins) + decoys
    oms.FASTAFile().store(output_path, combined)

    return len(decoys)


@click.command(help="Generate a target-decoy FASTA database.")
@click.option("--input", "input_path", required=True, help="Input target FASTA file")
@click.option("--output", "output_path", required=True, help="Output target-decoy FASTA file")
@click.option(
    "--method",
    type=click.Choice(["reverse", "shuffle"]),
    default="reverse",
    help="Decoy generation method (default: reverse)",
)
def main(input_path, output_path, method):
    count = generate_decoys(input_path, output_path, method=method)
    print(f"Generated {count} decoy proteins")
    print(f"Output written to {output_path}")


if __name__ == "__main__":
    main()
