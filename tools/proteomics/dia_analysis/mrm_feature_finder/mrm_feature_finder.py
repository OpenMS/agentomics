"""
MRM Feature Finder
==================
Find and score features in MRM/SRM chromatograms using MRMFeatureFinderScoring.

Features
--------
- Pick peaks in chromatographic data
- Score transition groups
- Output features as featureXML

Usage
-----
    python mrm_feature_finder.py --input chromatograms.mzML --transitions transitions.tsv --output features.featureXML
"""

import csv
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def load_transitions_to_targeted(transitions_path: str) -> oms.TargetedExperiment:
    """Load transitions TSV into a TargetedExperiment.

    Expected columns: PrecursorMz, ProductMz, LibraryIntensity,
    PeptideSequence, ProteinName, transition_name, transition_group_id.

    Parameters
    ----------
    transitions_path : str
        Path to transitions TSV file.

    Returns
    -------
    oms.TargetedExperiment
        Targeted experiment with transitions.
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


def find_mrm_features(
    input_path: str,
    transitions_path: str,
    output_path: str,
) -> int:
    """Find features in MRM chromatograms.

    Parameters
    ----------
    input_path : str
        Path to input chromatograms mzML file.
    transitions_path : str
        Path to transitions TSV file.
    output_path : str
        Path to output featureXML file.

    Returns
    -------
    int
        Number of features found.
    """
    chrom_exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, chrom_exp)

    te = load_transitions_to_targeted(transitions_path)

    mff = oms.MRMFeatureFinderScoring()
    params = mff.getParameters()
    params.setValue(b"TransitionGroupPicker:min_peak_width", 1.0)
    params.setValue(
        b"TransitionGroupPicker:PeakPickerChromatogram:signal_to_noise", 0.1
    )
    params.setValue(
        b"TransitionGroupPicker:PeakPickerChromatogram:gauss_width", 30.0
    )
    mff.setParameters(params)

    output = oms.FeatureMap()
    trafo = oms.TransformationDescription()
    swath_map = oms.MSExperiment()

    mff.pickExperiment(chrom_exp, output, te, trafo, swath_map)

    oms.FeatureXMLFile().store(output_path, output)
    return output.size()


@click.command(help="Find features in MRM/SRM chromatograms.")
@click.option("--input", "input_path", required=True, help="Path to chromatograms mzML.")
@click.option("--transitions", required=True, help="Path to transitions TSV file.")
@click.option("--output", required=True, help="Path to output featureXML file.")
def main(input_path, transitions, output):
    """CLI entry point."""
    n = find_mrm_features(input_path, transitions, output)
    print(f"Found {n} features -> {output}")


if __name__ == "__main__":
    main()
