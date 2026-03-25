"""
Mass Difference Network Builder
=================================
Connect features by known biotransformation mass differences to build
a molecular network.

Each edge represents a mass difference matching a known reaction
(e.g., oxidation, dehydration, methylation) within the given tolerance.

Usage
-----
    python mass_difference_network_builder.py --input features.tsv --reactions reactions.tsv \
        --tolerance 0.005 --output network.tsv
"""

import csv

import click
import pyopenms as oms  # noqa: F401

# Built-in common biotransformations: (name, mass_diff)
DEFAULT_REACTIONS = [
    ("Oxidation", 15.994915),
    ("Reduction", -15.994915),
    ("Dehydration", -18.010565),
    ("Hydration", 18.010565),
    ("Methylation", 14.015650),
    ("Demethylation", -14.015650),
    ("Acetylation", 42.010565),
    ("Deacetylation", -42.010565),
    ("Glucuronidation", 176.032088),
    ("Sulfation", 79.956815),
    ("Phosphorylation", 79.966331),
    ("Hydrogenation", 2.015650),
    ("Dehydrogenation", -2.015650),
    ("Hydroxylation", 15.994915),
    ("Glycine conjugation", 57.021464),
]


def load_reactions(path: str | None) -> list[tuple[str, float]]:
    """Load reaction definitions from a TSV file.

    Parameters
    ----------
    path:
        Path to TSV with columns: reaction_name, mass_diff.
        If None, returns the built-in default reactions.

    Returns
    -------
    list[tuple[str, float]]
    """
    if path is None:
        return DEFAULT_REACTIONS

    reactions = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            reactions.append((row["reaction_name"], float(row["mass_diff"])))
    return reactions


def build_network(
    features: list[dict],
    reactions: list[tuple[str, float]],
    tolerance: float = 0.005,
) -> list[dict]:
    """Build a mass-difference network from features.

    Parameters
    ----------
    features:
        List of dicts with at least keys: feature_id, mz.
    reactions:
        List of (reaction_name, mass_diff) tuples.
    tolerance:
        Absolute mass tolerance in Da.

    Returns
    -------
    list[dict]
        List of edges: source_id, target_id, reaction, mass_diff, error_da.
    """
    edges = []
    n = len(features)

    for i in range(n):
        mz_i = float(features[i]["mz"])
        id_i = features[i].get("feature_id", str(i))

        for j in range(i + 1, n):
            mz_j = float(features[j]["mz"])
            id_j = features[j].get("feature_id", str(j))
            diff = mz_j - mz_i

            for rxn_name, rxn_diff in reactions:
                error = abs(diff - rxn_diff)
                if error <= tolerance:
                    edges.append({
                        "source_id": id_i,
                        "target_id": id_j,
                        "reaction": rxn_name,
                        "mass_diff": round(diff, 6),
                        "error_da": round(error, 6),
                    })
                # Also check reverse direction
                error_rev = abs(-diff - rxn_diff)
                if error_rev <= tolerance:
                    edges.append({
                        "source_id": id_j,
                        "target_id": id_i,
                        "reaction": rxn_name,
                        "mass_diff": round(-diff, 6),
                        "error_da": round(error_rev, 6),
                    })

    return edges


@click.command()
@click.option("--input", "input_file", required=True, help="Features TSV (feature_id, mz)")
@click.option("--reactions", "reactions_file", default=None,
              help="Reactions TSV (reaction_name, mass_diff). Uses built-in table if omitted.")
@click.option("--tolerance", type=float, default=0.005,
              help="Mass tolerance in Da (default: 0.005)")
@click.option("--output", required=True, help="Output network TSV")
def main(input_file, reactions_file, tolerance, output):
    features = []
    with open(input_file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            features.append(row)

    reactions = load_reactions(reactions_file)
    edges = build_network(features, reactions, tolerance=tolerance)

    with open(output, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["source_id", "target_id", "reaction", "mass_diff", "error_da"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(edges)

    print(f"Network: {len(edges)} edges from {len(features)} features, written to {output}")


if __name__ == "__main__":
    main()
