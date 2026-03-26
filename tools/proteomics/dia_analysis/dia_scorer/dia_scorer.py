"""
DIA Scorer
==========
Score DIA spectra against transition lists using DIAScoring.

Features
--------
- Score DIA spectra for transition matches using isotope-based scoring
- Compute dot product and Manhattan distance scores
- Output scores as TSV

Usage
-----
    python dia_scorer.py --input dia.mzML --transitions transitions.tsv --output scores.tsv
"""

import csv
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def load_transitions(transitions_path: str) -> list:
    """Load transitions from TSV into a list of dicts.

    Parameters
    ----------
    transitions_path : str
        Path to transitions TSV file with columns PrecursorMz, ProductMz, etc.

    Returns
    -------
    list
        List of transition dicts.
    """
    transitions = []
    with open(transitions_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            transitions.append({
                "PrecursorMz": float(row["PrecursorMz"]),
                "ProductMz": float(row["ProductMz"]),
                "LibraryIntensity": float(row.get("LibraryIntensity", 100.0)),
                "PeptideSequence": row.get("PeptideSequence", ""),
                "ProteinName": row.get("ProteinName", ""),
                "transition_name": row.get("transition_name", ""),
                "transition_group_id": row.get(
                    "transition_group_id", row.get("PeptideSequence", "")
                ),
            })
    return transitions


def _group_transitions(transitions: list) -> dict:
    """Group transitions by transition_group_id.

    Parameters
    ----------
    transitions : list
        List of transition dicts.

    Returns
    -------
    dict
        Mapping from group_id to list of transitions.
    """
    groups = {}
    for t in transitions:
        gid = t["transition_group_id"]
        if gid not in groups:
            groups[gid] = []
        groups[gid].append(t)
    return groups


def score_dia(
    input_path: str,
    transitions_path: str,
    output_path: str,
) -> int:
    """Score DIA spectra against transitions using DIAScoring.

    For each transition group, finds the best-matching spectrum and computes
    isotope-based dot product and Manhattan distance scores.

    Parameters
    ----------
    input_path : str
        Path to input DIA mzML file.
    transitions_path : str
        Path to transitions TSV file.
    output_path : str
        Path to output scores TSV file.

    Returns
    -------
    int
        Number of scored transition groups.
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    transitions = load_transitions(transitions_path)
    groups = _group_transitions(transitions)

    scorer = oms.DIAScoring()
    params = scorer.getParameters()
    params.setValue(b"dia_extraction_window", 0.05)
    params.setValue(b"dia_centroided", "false")
    scorer.setParameters(params)

    results = []
    spectra = []
    for i in range(exp.getNrSpectra()):
        spectra.append(exp.getSpectrum(i))

    for group_id, group_trans in groups.items():
        light_transitions = []
        for t in group_trans:
            lt = oms.LightTransition()
            lt.precursor_mz = t["PrecursorMz"]
            lt.product_mz = t["ProductMz"]
            lt.library_intensity = t["LibraryIntensity"]
            lt.transition_name = t["transition_name"].encode()
            lt.peptide_ref = group_id.encode()
            lt.detecting_transition = True
            lt.quantifying_transition = True
            lt.identifying_transition = False
            light_transitions.append(lt)

        best_dotprod = 0.0
        best_manhattan = 1.0
        best_rt = 0.0
        scored = False

        for spec in spectra:
            if spec.getMSLevel() != 2:
                continue

            os_spectra = [oms.OSSpectrum()]
            mzs_arr, ints_arr = spec.get_peaks()
            if len(mzs_arr) == 0:
                continue
            os_spectra[0].set_mz_array(list(mzs_arr))
            os_spectra[0].set_intensity_array(list(ints_arr))

            im_range = oms.RangeMobility()

            try:
                dotprod = [0.0]
                manhattan = [0.0]
                scorer.score_with_isotopes(
                    os_spectra, light_transitions, im_range,
                    dotprod, manhattan,
                )
                if dotprod[0] > best_dotprod:
                    best_dotprod = dotprod[0]
                    best_manhattan = manhattan[0]
                    best_rt = spec.getRT()
                    scored = True
            except (RuntimeError, Exception):
                # Some spectra may not score well; continue
                pass

        results.append({
            "transition_group_id": group_id,
            "peptide_sequence": group_trans[0]["PeptideSequence"],
            "protein_name": group_trans[0]["ProteinName"],
            "precursor_mz": group_trans[0]["PrecursorMz"],
            "n_transitions": len(group_trans),
            "best_rt": round(best_rt, 4),
            "dotprod_score": round(best_dotprod, 6),
            "manhattan_score": round(best_manhattan, 6),
            "scored": scored,
        })

    fieldnames = [
        "transition_group_id", "peptide_sequence", "protein_name",
        "precursor_mz", "n_transitions", "best_rt",
        "dotprod_score", "manhattan_score", "scored",
    ]
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    return len(results)


@click.command(help="Score DIA spectra against a transition list.")
@click.option("--input", "input_path", required=True, help="Path to input DIA mzML file.")
@click.option("--transitions", required=True, help="Path to transitions TSV file.")
@click.option("--output", required=True, help="Path to output scores TSV file.")
def main(input_path, transitions, output):
    """CLI entry point."""
    n = score_dia(input_path, transitions, output)
    print(f"Scored {n} transition groups -> {output}")


if __name__ == "__main__":
    main()
