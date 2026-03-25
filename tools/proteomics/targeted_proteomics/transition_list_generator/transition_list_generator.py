"""
Transition List Generator
==========================
Generate SRM/MRM/PRM transition lists from peptide sequences.

Features
--------
- Generate b and y fragment ion transitions
- Support multiple charge states
- Filter product ions by ion series and range
- Output in standard transition list format

Usage
-----
    python transition_list_generator.py --peptides PEPTIDEK,ANOTHERPEPTIDE --charge 2,3 --output transitions.tsv
    python transition_list_generator.py --peptides PEPTIDEK --charge 2 --product-ions y3-y8 --output transitions.tsv
"""

import csv

import click
import pyopenms as oms

PROTON = 1.007276


def parse_ion_range(ion_range: str) -> tuple:
    """Parse an ion range string like 'y3-y8' into (series, start, end).

    Parameters
    ----------
    ion_range : str
        Ion range string (e.g., 'y3-y8', 'b2-b6').

    Returns
    -------
    tuple
        (series_letter, start_number, end_number).
    """
    parts = ion_range.split("-")
    series = parts[0][0]
    start = int(parts[0][1:])
    end = int(parts[1][1:]) if len(parts) > 1 else start
    return series, start, end


def generate_transitions(sequence: str, precursor_charges: list = None,
                         product_charges: list = None,
                         ion_range: str = None) -> list:
    """Generate SRM/MRM transitions for a peptide.

    Parameters
    ----------
    sequence : str
        Peptide sequence.
    precursor_charges : list
        List of precursor charge states.
    product_charges : list
        List of product ion charge states.
    ion_range : str
        Optional filter like 'y3-y8' to limit product ions.

    Returns
    -------
    list
        List of transition dicts.
    """
    if precursor_charges is None:
        precursor_charges = [2]
    if product_charges is None:
        product_charges = [1]

    aa_seq = oms.AASequence.fromString(sequence)
    mono_mass = aa_seq.getMonoWeight()

    # Generate theoretical spectrum
    tsg = oms.TheoreticalSpectrumGenerator()
    params = tsg.getParameters()
    params.setValue("add_b_ions", "true")
    params.setValue("add_y_ions", "true")
    params.setValue("add_metainfo", "true")
    tsg.setParameters(params)

    transitions = []
    for prec_z in precursor_charges:
        prec_mz = (mono_mass + prec_z * PROTON) / prec_z

        for prod_z in product_charges:
            spec = oms.MSSpectrum()
            tsg.getSpectrum(spec, aa_seq, prod_z, prec_z)

            for i in range(spec.size()):
                peak = spec[i]
                product_mz = peak.getMZ()
                annotation = ""
                if spec.getStringDataArrays():
                    ann = spec.getStringDataArrays()[0][i]
                    annotation = ann.decode() if isinstance(ann, bytes) else str(ann)

                # Apply ion range filter if specified
                if ion_range:
                    series, start, end = parse_ion_range(ion_range)
                    if not annotation.startswith(series):
                        continue
                    # Extract ion number from annotation
                    ion_num_str = ""
                    for ch in annotation[1:]:
                        if ch.isdigit():
                            ion_num_str += ch
                        else:
                            break
                    if ion_num_str:
                        ion_num = int(ion_num_str)
                        if ion_num < start or ion_num > end:
                            continue

                transitions.append({
                    "peptide": sequence,
                    "precursor_mz": round(prec_mz, 6),
                    "precursor_charge": prec_z,
                    "product_mz": round(product_mz, 6),
                    "product_charge": prod_z,
                    "annotation": annotation,
                })

    return transitions


@click.command(help="Generate SRM/MRM/PRM transition lists.")
@click.option("--peptides", required=True, help="Comma-separated peptide sequences.")
@click.option("--charge", type=str, default="2", help="Comma-separated precursor charges (default: 2).")
@click.option("--product-charge", type=str, default="1", help="Comma-separated product ion charges (default: 1).")
@click.option("--product-ions", type=str, default=None, help="Ion range filter (e.g., 'y3-y8').")
@click.option("--output", default=None, help="Output TSV file.")
def main(peptides, charge, product_charge, product_ions, output):
    """CLI entry point."""
    peptide_list = [p.strip() for p in peptides.split(",") if p.strip()]
    precursor_charges = [int(c.strip()) for c in charge.split(",")]
    product_charges_list = [int(c.strip()) for c in product_charge.split(",")]

    all_transitions = []
    for pep in peptide_list:
        transitions = generate_transitions(pep, precursor_charges, product_charges_list, product_ions)
        all_transitions.extend(transitions)

    if output:
        with open(output, "w", newline="") as fh:
            fieldnames = ["peptide", "precursor_mz", "precursor_charge", "product_mz",
                          "product_charge", "annotation"]
            writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(all_transitions)
        print(f"Generated {len(all_transitions)} transitions -> {output}")
    else:
        for t in all_transitions:
            print(f"{t['peptide']}\t{t['precursor_mz']}\t{t['product_mz']}\t{t['annotation']}")


if __name__ == "__main__":
    main()
