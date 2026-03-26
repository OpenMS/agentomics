"""
MRM Transition Group Picker
============================
Pick peaks in MRM transition group chromatograms using MRMTransitionGroupPicker.

Features
--------
- Group chromatograms by transition group
- Pick peaks across co-eluting transitions
- Output picked chromatograms with feature annotations

Usage
-----
    python mrm_transition_group_picker.py --input chromatograms.mzML --output picked.mzML
    python mrm_transition_group_picker.py --input chromatograms.mzML --output picked.mzML --min-peak-width 5.0
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def _parse_group_id(native_id: str) -> str:
    """Extract transition group ID from a chromatogram native ID.

    Convention: native IDs of the form 'groupid_transitionindex' or just
    the native ID itself as the group.

    Parameters
    ----------
    native_id : str
        Chromatogram native ID.

    Returns
    -------
    str
        Transition group ID.
    """
    # If the native ID has the form 'groupid_N', group by prefix
    parts = native_id.rsplit("_", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return native_id


def pick_transition_groups(
    input_path: str,
    output_path: str,
    min_peak_width: float = 0.001,
    signal_to_noise: float = 0.1,
) -> int:
    """Pick peaks in transition group chromatograms.

    Groups chromatograms by transition group (based on native ID convention)
    and uses MRMTransitionGroupPicker to find consensus peaks.

    Parameters
    ----------
    input_path : str
        Path to input chromatograms mzML file.
    output_path : str
        Path to output mzML file with picked features.
    min_peak_width : float
        Minimum peak width in RT units.
    signal_to_noise : float
        Signal-to-noise threshold.

    Returns
    -------
    int
        Number of features picked across all transition groups.
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    tgp = oms.MRMTransitionGroupPicker()
    params = tgp.getParameters()
    params.setValue(b"min_peak_width", min_peak_width)
    params.setValue(b"PeakPickerChromatogram:signal_to_noise", signal_to_noise)
    tgp.setParameters(params)

    # Group chromatograms by transition group
    groups = {}
    for i in range(exp.getNrChromatograms()):
        chrom = exp.getChromatogram(i)
        nid = chrom.getNativeID()
        if isinstance(nid, bytes):
            nid = nid.decode()
        group_id = _parse_group_id(nid)
        if group_id not in groups:
            groups[group_id] = []
        groups[group_id].append(chrom)

    total_features = 0
    out_exp = oms.MSExperiment()

    for group_id, chroms in groups.items():
        tg = oms.MRMTransitionGroupCP()
        tg.setTransitionGroupID(group_id.encode())

        for chrom in chroms:
            nid = chrom.getNativeID()
            if isinstance(nid, bytes):
                nid_str = nid
            else:
                nid_str = nid.encode()

            t = oms.ReactionMonitoringTransition()
            t.setNativeID(nid_str)
            t.setDetectingTransition(True)
            t.setQuantifyingTransition(True)
            t.setIdentifyingTransition(False)
            tg.addTransition(t, nid_str)
            tg.addChromatogram(chrom, nid_str)

        tgp.pickTransitionGroup(tg)
        features = tg.getFeatures()
        total_features += len(features)

        # Store the chromatograms (unchanged) to output
        for chrom in chroms:
            out_exp.addChromatogram(chrom)

    oms.MzMLFile().store(output_path, out_exp)
    return total_features


@click.command(help="Pick peaks in MRM transition group chromatograms.")
@click.option("--input", "input_path", required=True, help="Path to chromatograms mzML.")
@click.option("--output", required=True, help="Path to output mzML file.")
@click.option("--min-peak-width", type=float, default=0.001, help="Minimum peak width.")
@click.option("--signal-to-noise", type=float, default=0.1, help="S/N threshold.")
def main(input_path, output, min_peak_width, signal_to_noise):
    """CLI entry point."""
    n = pick_transition_groups(input_path, output, min_peak_width, signal_to_noise)
    print(f"Picked {n} features -> {output}")


if __name__ == "__main__":
    main()
