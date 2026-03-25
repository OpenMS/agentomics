"""
Spectral Library Format Converter
==================================
Convert between spectral library formats (MSP to TraML).
Parses MSP format and creates a TraML TargetedExperiment.

Features:
- Parse MSP spectral library files
- Convert to TraML format
- Preserve peptide metadata and transitions

Usage
-----
    python spectral_library_format_converter.py --input library.msp --output library.traml --format traml
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def parse_msp(filepath: str) -> list[dict]:
    """Parse an MSP spectral library file.

    Parameters
    ----------
    filepath : str
        Path to MSP file.

    Returns
    -------
    list[dict]
        List of spectrum dicts with keys: name, mw, comment, precursor_mz,
        charge, peaks (list of (mz, intensity) tuples).
    """
    spectra = []
    current = None

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                if current is not None and current.get("peaks"):
                    spectra.append(current)
                    current = None
                continue

            if line.startswith("Name:"):
                if current is not None and current.get("peaks"):
                    spectra.append(current)
                current = {
                    "name": line[5:].strip(),
                    "mw": 0.0,
                    "comment": "",
                    "precursor_mz": 0.0,
                    "charge": 0,
                    "peaks": [],
                }
            elif current is not None:
                if line.startswith("MW:"):
                    current["mw"] = float(line[3:].strip())
                elif line.startswith("Comment:"):
                    current["comment"] = line[8:].strip()
                    # Extract charge from comment if present
                    for part in current["comment"].split():
                        if part.startswith("Charge="):
                            current["charge"] = int(part.split("=")[1].replace("+", "").replace("-", ""))
                        elif part.startswith("Parent="):
                            current["precursor_mz"] = float(part.split("=")[1])
                elif line.startswith("Num peaks:") or line.startswith("Num Peaks:"):
                    pass  # Skip, we count peaks from data
                elif line[0].isdigit():
                    parts = line.split()
                    if len(parts) >= 2:
                        mz = float(parts[0])
                        intensity = float(parts[1])
                        current["peaks"].append((mz, intensity))

    if current is not None and current.get("peaks"):
        spectra.append(current)

    return spectra


def msp_to_targeted_experiment(spectra: list[dict]) -> oms.TargetedExperiment:
    """Convert parsed MSP spectra to a TargetedExperiment.

    Parameters
    ----------
    spectra : list[dict]
        List of parsed MSP spectrum dictionaries.

    Returns
    -------
    oms.TargetedExperiment
        pyopenms TargetedExperiment object.
    """
    targeted_exp = oms.TargetedExperiment()
    proteins = []
    peptides = []
    transitions = []

    for spec_idx, spec in enumerate(spectra):
        peptide_id = f"peptide_{spec_idx}".encode()
        protein_id = f"protein_{spec_idx}".encode()

        # Create protein
        protein = oms.Protein()
        protein.id = protein_id
        proteins.append(protein)

        # Create peptide
        peptide = oms.Peptide()
        peptide.id = peptide_id
        peptide.protein_refs = [protein_id]
        seq_str = spec["name"].split("/")[0] if "/" in spec["name"] else spec["name"]
        peptide.sequence = seq_str.encode() if isinstance(seq_str, str) else seq_str
        if spec["charge"] > 0:
            peptide.setChargeState(spec["charge"])
        peptides.append(peptide)

        # Create transitions from peaks
        precursor_mz = spec["precursor_mz"]
        if precursor_mz == 0 and spec["mw"] > 0 and spec["charge"] > 0:
            precursor_mz = (spec["mw"] + spec["charge"] * 1.007276) / spec["charge"]

        for peak_idx, (mz, intensity) in enumerate(spec["peaks"]):
            transition = oms.ReactionMonitoringTransition()
            transition.setNativeID(f"transition_{spec_idx}_{peak_idx}".encode())
            transition.setPeptideRef(peptide_id)
            transition.setPrecursorMZ(precursor_mz)
            transition.setProductMZ(mz)
            transition.setLibraryIntensity(intensity)
            transitions.append(transition)

    targeted_exp.setProteins(proteins)
    targeted_exp.setPeptides(peptides)
    targeted_exp.setTransitions(transitions)

    return targeted_exp


def convert_msp_to_traml(input_path: str, output_path: str) -> int:
    """Convert MSP file to TraML format.

    Parameters
    ----------
    input_path : str
        Path to input MSP file.
    output_path : str
        Path to output TraML file.

    Returns
    -------
    int
        Number of spectra converted.
    """
    spectra = parse_msp(input_path)
    targeted_exp = msp_to_targeted_experiment(spectra)
    oms.TraMLFile().store(output_path, targeted_exp)
    return len(spectra)


def create_synthetic_msp(output_path: str, n_spectra: int = 3) -> None:
    """Create a synthetic MSP file for testing.

    Parameters
    ----------
    output_path : str
        Path to write the synthetic MSP file.
    n_spectra : int
        Number of spectra to generate.
    """
    sequences = ["PEPTIDEK", "ACDEFGHIK", "LMNPQR"]

    with open(output_path, "w") as f:
        for i in range(n_spectra):
            seq = sequences[i % len(sequences)]
            charge = 2
            mw = 500.0 + i * 100
            precursor_mz = (mw + charge * 1.007276) / charge

            f.write(f"Name: {seq}/{charge}\n")
            f.write(f"MW: {mw}\n")
            f.write(f"Comment: Charge={charge}+ Parent={precursor_mz:.4f}\n")
            f.write("Num peaks: 5\n")
            for j in range(5):
                frag_mz = 100.0 + j * 50 + i * 10
                intensity = 1000.0 - j * 150
                f.write(f"{frag_mz:.4f}\t{intensity:.1f}\n")
            f.write("\n")


@click.command(help="Convert between spectral library formats (MSP to TraML).")
@click.option("--input", "input", required=True, help="Path to input MSP file")
@click.option("--output", required=True, help="Path to output file")
@click.option(
    "--format", "format", default="traml",
    type=click.Choice(["traml"]), help="Output format (default: traml)",
)
def main(input, output, format):
    count = convert_msp_to_traml(input, output)
    print(f"Converted {count} spectra to {format} format: {output}")


if __name__ == "__main__":
    main()
