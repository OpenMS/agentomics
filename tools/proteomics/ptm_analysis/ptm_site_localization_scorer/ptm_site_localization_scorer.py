"""
PTM Site Localization Scorer
==============================
Score PTM site localization confidence using fragment ion coverage comparison.

Features
--------
- Generate theoretical spectra for candidate PTM site assignments
- Match experimental peaks against theoretical fragments
- Score each candidate site by fragment ion coverage
- Report site localization probabilities

Usage
-----
    python ptm_site_localization_scorer.py --mz-list "200.1,300.2,400.3" --intensities "100,200,150" \\
        --peptide "PEPS(Phospho)TIDEK" --tolerance 0.02 --output scores.tsv
"""

import csv
import json

import click
import pyopenms as oms


def generate_theoretical_spectrum(sequence: str, charge: int = 1) -> list:
    """Generate theoretical b/y ion spectrum for a modified peptide.

    Parameters
    ----------
    sequence : str
        Modified peptide sequence in pyopenms notation.
    charge : int
        Precursor charge state.

    Returns
    -------
    list
        List of (mz, annotation) tuples.
    """
    aa_seq = oms.AASequence.fromString(sequence)
    spec = oms.MSSpectrum()
    tsg = oms.TheoreticalSpectrumGenerator()

    params = tsg.getParameters()
    params.setValue("add_b_ions", "true")
    params.setValue("add_y_ions", "true")
    params.setValue("add_metainfo", "true")
    tsg.setParameters(params)

    tsg.getSpectrum(spec, aa_seq, 1, charge)

    ions = []
    for i in range(spec.size()):
        peak = spec[i]
        mz = peak.getMZ()
        name = ""
        if spec.getStringDataArrays():
            name = spec.getStringDataArrays()[0][i].decode() if isinstance(
                spec.getStringDataArrays()[0][i], bytes
            ) else str(spec.getStringDataArrays()[0][i])
        ions.append((mz, name))
    return ions


def match_peaks(experimental_mz: list, theoretical_ions: list,
                tolerance: float = 0.02) -> list:
    """Match experimental peaks to theoretical fragment ions.

    Parameters
    ----------
    experimental_mz : list
        List of experimental m/z values.
    theoretical_ions : list
        List of (mz, annotation) tuples from theoretical spectrum.
    tolerance : float
        Mass tolerance in Da for peak matching.

    Returns
    -------
    list
        List of matched (experimental_mz, theoretical_mz, annotation, error) tuples.
    """
    matches = []
    for exp_mz in experimental_mz:
        for theo_mz, annotation in theoretical_ions:
            error = abs(exp_mz - theo_mz)
            if error <= tolerance:
                matches.append((exp_mz, theo_mz, annotation, round(error, 6)))
    return matches


def generate_site_candidates(sequence: str, mod_name: str, applicable_residues: str) -> list:
    """Generate all possible site assignment candidates for a modification.

    Parameters
    ----------
    sequence : str
        Unmodified peptide sequence.
    mod_name : str
        Modification name (e.g., 'Phospho').
    applicable_residues : str
        String of residues that can carry the modification (e.g., 'STY').

    Returns
    -------
    list
        List of modified sequence strings, one per candidate site.
    """
    plain = oms.AASequence.fromString(sequence).toUnmodifiedString()
    candidates = []
    for i, aa in enumerate(plain):
        if aa in applicable_residues:
            seq_list = list(plain)
            seq_list[i] = f"{aa}({mod_name})"
            candidates.append("".join(seq_list))
    return candidates


def score_localization(experimental_mz: list, experimental_intensities: list,
                       peptide: str, tolerance: float = 0.02,
                       charge: int = 1) -> dict:
    """Score PTM site localization for a modified peptide.

    Parameters
    ----------
    experimental_mz : list
        Experimental m/z values.
    experimental_intensities : list
        Experimental intensities corresponding to m/z values.
    peptide : str
        Modified peptide sequence.
    tolerance : float
        Mass tolerance for matching.
    charge : int
        Charge state.

    Returns
    -------
    dict
        Dictionary with the input peptide score and candidate site scores.
    """
    # Score the given assignment
    theo_ions = generate_theoretical_spectrum(peptide, charge)
    matches = match_peaks(experimental_mz, theo_ions, tolerance)
    given_score = len(matches)

    # Generate alternative candidates
    aa_seq = oms.AASequence.fromString(peptide)
    plain = aa_seq.toUnmodifiedString()

    # Find the modification in the given peptide
    mod_info = None
    for i in range(aa_seq.size()):
        if aa_seq.getResidue(i).isModified():
            mod_info = (i, aa_seq.getResidue(i).getModificationName())
            break

    candidate_scores = []
    if mod_info:
        mod_pos, mod_name = mod_info
        # Determine applicable residues from the mod origin
        applicable = set()
        mod_db = oms.ModificationsDB()
        mods_set = set()
        mod_db.searchModifications(mods_set, mod_name, "", 0)
        for mod_obj in mods_set:
            origin = mod_obj.getOrigin()
            if origin:
                applicable.add(origin)

        if not applicable:
            applicable = {plain[mod_pos]}

        candidates = generate_site_candidates(plain, mod_name, "".join(applicable))
        total_matches = 0
        for cand_seq in candidates:
            cand_ions = generate_theoretical_spectrum(cand_seq, charge)
            cand_matches = match_peaks(experimental_mz, cand_ions, tolerance)
            cand_score = len(cand_matches)
            total_matches += max(cand_score, 1)
            candidate_scores.append({
                "sequence": cand_seq,
                "matched_ions": cand_score,
            })

        # Normalize to probabilities
        for cs in candidate_scores:
            cs["probability"] = round(cs["matched_ions"] / total_matches, 4) if total_matches > 0 else 0.0
    else:
        candidate_scores.append({
            "sequence": peptide,
            "matched_ions": given_score,
            "probability": 1.0,
        })

    candidate_scores.sort(key=lambda x: x["matched_ions"], reverse=True)

    return {
        "peptide": peptide,
        "charge": charge,
        "tolerance": tolerance,
        "total_experimental_peaks": len(experimental_mz),
        "matched_ions_given": given_score,
        "candidates": candidate_scores,
    }


@click.command(help="Score PTM site localization confidence.")
@click.option("--mz-list", required=True, help="Comma-separated experimental m/z values.")
@click.option("--intensities", required=True, help="Comma-separated intensities.")
@click.option("--peptide", required=True, help="Modified peptide sequence.")
@click.option("--tolerance", type=float, default=0.02, help="Mass tolerance in Da (default: 0.02).")
@click.option("--charge", type=int, default=1, help="Charge state (default: 1).")
@click.option("--output", default=None, help="Output file (.tsv or .json).")
def main(mz_list, intensities, peptide, tolerance, charge, output):
    """CLI entry point."""
    mz_values = [float(x.strip()) for x in mz_list.split(",")]
    intensities_list = [float(x.strip()) for x in intensities.split(",")]

    result = score_localization(mz_values, intensities_list, peptide, tolerance, charge)

    if output:
        if output.endswith(".json"):
            with open(output, "w") as fh:
                json.dump(result, fh, indent=2)
        else:
            with open(output, "w", newline="") as fh:
                writer = csv.DictWriter(
                    fh, fieldnames=["sequence", "matched_ions", "probability"], delimiter="\t"
                )
                writer.writeheader()
                writer.writerows(result["candidates"])
        print(f"Results written to {output}")
    else:
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
