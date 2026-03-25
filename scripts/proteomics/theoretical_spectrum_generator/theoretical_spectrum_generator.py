"""
Theoretical Spectrum Generator
==============================
Generate theoretical b/y/a/c/x/z fragment ion spectra for a peptide sequence.
Output annotated TSV with ion type, ion number, charge, m/z, and annotation.

Features:
- Supports b, y, a, c, x, z ion types
- Multiple charge states
- Optional neutral losses (H2O, NH3)
- TSV output with full annotation

Usage
-----
    python theoretical_spectrum_generator.py --sequence PEPTIDEK --charge 2 --ion-types b,y
    python theoretical_spectrum_generator.py --sequence PEPTIDEK --charge 1 \
        --ion-types b,y,a --add-losses --output fragments.tsv
"""

import argparse
import csv
import sys

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def generate_theoretical_spectrum(
    sequence: str,
    charge: int = 1,
    ion_types: list[str] | None = None,
    add_losses: bool = False,
) -> list[dict]:
    """Generate theoretical fragment ion spectrum for a peptide sequence.

    Parameters
    ----------
    sequence : str
        Amino acid sequence, e.g. ``"PEPTIDEK"``.
    charge : int
        Maximum charge state for fragment ions (default 1).
    ion_types : list[str] or None
        Ion types to generate, e.g. ``["b", "y"]``. Defaults to ``["b", "y"]``.
    add_losses : bool
        Whether to include neutral losses (H2O, NH3).

    Returns
    -------
    list[dict]
        List of dicts with keys: ion_type, ion_number, charge, mz, annotation.
    """
    if ion_types is None:
        ion_types = ["b", "y"]

    aa_seq = oms.AASequence.fromString(sequence)
    spec = oms.MSSpectrum()
    tsg = oms.TheoreticalSpectrumGenerator()

    param = tsg.getParameters()
    param.setValue("add_b_ions", "true" if "b" in ion_types else "false")
    param.setValue("add_y_ions", "true" if "y" in ion_types else "false")
    param.setValue("add_a_ions", "true" if "a" in ion_types else "false")
    param.setValue("add_c_ions", "true" if "c" in ion_types else "false")
    param.setValue("add_x_ions", "true" if "x" in ion_types else "false")
    param.setValue("add_z_ions", "true" if "z" in ion_types else "false")
    param.setValue("add_losses", "true" if add_losses else "false")
    param.setValue("add_metainfo", "true")
    tsg.setParameters(param)

    tsg.getSpectrum(spec, aa_seq, charge, charge)

    results = []
    for i in range(spec.size()):
        mz = spec[i].getMZ()
        annotation_str = ""
        ion_type_str = ""
        ion_number = 0
        ion_charge = 1

        if spec.getStringDataArrays():
            annotations = spec.getStringDataArrays()
            if len(annotations) > 0 and i < annotations[0].size():
                annotation_str = annotations[0][i].decode() if isinstance(annotations[0][i], bytes) else str(
                    annotations[0][i]
                )

        if annotation_str:
            parsed = _parse_annotation(annotation_str)
            ion_type_str = parsed["ion_type"]
            ion_number = parsed["ion_number"]
            ion_charge = parsed["charge"]
        else:
            ion_type_str = "unknown"

        results.append({
            "ion_type": ion_type_str,
            "ion_number": ion_number,
            "charge": ion_charge,
            "mz": round(mz, 6),
            "annotation": annotation_str if annotation_str else f"{mz:.4f}",
        })

    return results


def _parse_annotation(annotation: str) -> dict:
    """Parse a pyopenms ion annotation string.

    Parameters
    ----------
    annotation : str
        Annotation string like ``"y3++"`` or ``"b5-H2O+"``.

    Returns
    -------
    dict
        Dictionary with ion_type, ion_number, charge.
    """
    ion_type = ""
    ion_number = 0
    charge = 1

    if not annotation:
        return {"ion_type": "unknown", "ion_number": 0, "charge": 1}

    # Extract ion type (first letter)
    for c in annotation:
        if c.isalpha():
            ion_type = c
            break

    # Extract ion number (digits after ion type letter)
    num_str = ""
    started = False
    for c in annotation:
        if c.isdigit():
            num_str += c
            started = True
        elif started:
            break
    if num_str:
        ion_number = int(num_str)

    # Count charge from + signs
    charge = annotation.count("+")
    if charge == 0:
        charge = 1

    return {"ion_type": ion_type, "ion_number": ion_number, "charge": charge}


def write_tsv(results: list[dict], output_path: str) -> None:
    """Write spectrum results to TSV file.

    Parameters
    ----------
    results : list[dict]
        List of fragment ion dictionaries.
    output_path : str
        Path to output TSV file.
    """
    fieldnames = ["ion_type", "ion_number", "charge", "mz", "annotation"]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


def main():
    parser = argparse.ArgumentParser(
        description="Generate theoretical fragment ion spectra for a peptide sequence."
    )
    parser.add_argument("--sequence", required=True, help="Amino acid sequence (e.g. PEPTIDEK)")
    parser.add_argument("--charge", type=int, default=1, help="Max charge state for fragment ions (default: 1)")
    parser.add_argument("--ion-types", default="b,y", help="Comma-separated ion types: b,y,a,c,x,z (default: b,y)")
    parser.add_argument("--add-losses", action="store_true", help="Include neutral losses (H2O, NH3)")
    parser.add_argument("--output", default=None, help="Output TSV file path (default: print to stdout)")
    args = parser.parse_args()

    ion_types = [t.strip() for t in args.ion_types.split(",")]
    results = generate_theoretical_spectrum(args.sequence, args.charge, ion_types, args.add_losses)

    if args.output:
        write_tsv(results, args.output)
        print(f"Wrote {len(results)} fragment ions to {args.output}")
    else:
        print(f"{'ion_type'}\t{'ion_number'}\t{'charge'}\t{'mz'}\t{'annotation'}")
        for r in results:
            print(f"{r['ion_type']}\t{r['ion_number']}\t{r['charge']}\t{r['mz']}\t{r['annotation']}")


if __name__ == "__main__":
    main()
