"""
Chromatogram Extractor
======================
Extract ion chromatograms from DIA/SWATH-MS mzML files using transition lists.

Features
--------
- Extract XICs for target transitions from DIA data
- Configurable m/z tolerance (ppm)
- Outputs chromatograms as mzML

Usage
-----
    python chromatogram_extractor.py --input dia.mzML --transitions transitions.tsv --output chromatograms.mzML
    python chromatogram_extractor.py --input dia.mzML --transitions transitions.tsv \
        --output chromatograms.mzML --mz-tol 10.0
"""

import csv
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def load_transitions_tsv(transitions_path: str) -> oms.TargetedExperiment:
    """Load transitions from a TSV file into a TargetedExperiment.

    Expected TSV columns: PrecursorMz, ProductMz, LibraryIntensity,
    PeptideSequence, ProteinName, transition_name, transition_group_id.

    Parameters
    ----------
    transitions_path : str
        Path to transitions TSV file.

    Returns
    -------
    oms.TargetedExperiment
        Targeted experiment with transitions loaded.
    """
    te = oms.TargetedExperiment()
    proteins_seen = set()
    peptides_seen = set()

    with open(transitions_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            protein_name = row.get("ProteinName", "UnknownProtein")
            if protein_name not in proteins_seen:
                prot = oms.Protein()
                prot.id = protein_name.encode()
                te.addProtein(prot)
                proteins_seen.add(protein_name)

            group_id = row.get("transition_group_id", row.get("PeptideSequence", ""))
            peptide_seq = row.get("PeptideSequence", "")
            if group_id not in peptides_seen:
                pep = oms.Peptide()
                pep.id = group_id.encode()
                pep.sequence = peptide_seq.encode()
                pep.protein_refs = [protein_name.encode()]
                charge = int(row.get("PrecursorCharge", 2))
                pep.setChargeState(charge)
                te.addPeptide(pep)
                peptides_seen.add(group_id)

            t = oms.ReactionMonitoringTransition()
            native_id = row.get("transition_name", f"tr_{len(te.getTransitions())}")
            t.setNativeID(native_id.encode())
            t.setPeptideRef(group_id.encode())
            t.setPrecursorMZ(float(row["PrecursorMz"]))
            t.setProductMZ(float(row["ProductMz"]))
            t.setLibraryIntensity(float(row.get("LibraryIntensity", 100.0)))
            t.setDetectingTransition(True)
            t.setQuantifyingTransition(True)
            t.setIdentifyingTransition(False)
            te.addTransition(t)

    return te


def extract_chromatograms(
    input_path: str,
    transitions_path: str,
    output_path: str,
    mz_tol: float = 10.0,
) -> int:
    """Extract chromatograms from a DIA mzML using transition definitions.

    Parameters
    ----------
    input_path : str
        Path to input DIA mzML file.
    transitions_path : str
        Path to transitions TSV file.
    output_path : str
        Path to output mzML file with chromatograms.
    mz_tol : float
        Mass tolerance in ppm (default 10.0).

    Returns
    -------
    int
        Number of chromatograms extracted.
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    te = load_transitions_tsv(transitions_path)

    ce = oms.ChromatogramExtractor()
    output_chroms = []
    coords = []
    ce.prepare_coordinates(output_chroms, coords, te, -1.0, False, 0)

    sao = oms.SpectrumAccessOpenMS(exp)
    ce.extractChromatograms(sao, output_chroms, coords, mz_tol, True, -1.0, b"tophat")

    out_exp = oms.MSExperiment()
    transitions = te.getTransitions()
    for i, oc in enumerate(output_chroms):
        chrom = oms.MSChromatogram()
        times = oc.get_time_array()
        intensities = oc.get_intensity_array()
        if len(times) > 0:
            chrom.set_peaks((list(times), list(intensities)))
        if i < len(transitions):
            chrom.setNativeID(transitions[i].getNativeID())
        else:
            chrom.setNativeID(f"chrom_{i}".encode())
        out_exp.addChromatogram(chrom)

    oms.MzMLFile().store(output_path, out_exp)
    return len(output_chroms)


@click.command(help="Extract chromatograms from DIA mzML using transition list.")
@click.option("--input", "input_path", required=True, help="Path to input DIA mzML file.")
@click.option("--transitions", required=True, help="Path to transitions TSV file.")
@click.option("--output", required=True, help="Path to output chromatograms mzML file.")
@click.option("--mz-tol", type=float, default=10.0, help="Mass tolerance in ppm (default: 10.0).")
def main(input_path, transitions, output, mz_tol):
    """CLI entry point."""
    n = extract_chromatograms(input_path, transitions, output, mz_tol)
    print(f"Extracted {n} chromatograms -> {output}")


if __name__ == "__main__":
    main()
