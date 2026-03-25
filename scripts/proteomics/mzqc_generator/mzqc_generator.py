"""
mzQC Generator
==============
Generate an mzQC-format JSON file from an mzML run. The output follows
a simplified mzQC schema with standard metric vocabulary identifiers.

Usage
-----
    python mzqc_generator.py --input run.mzML --output qc.mzQC
"""

import argparse
import json
import math
import sys
from datetime import datetime, timezone

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def generate_mzqc(exp: oms.MSExperiment, input_file: str = "unknown.mzML") -> dict:
    """Build an mzQC-style dictionary from an MSExperiment.

    Parameters
    ----------
    exp:
        Loaded ``pyopenms.MSExperiment`` instance.
    input_file:
        Name of the source mzML file (used in metadata).

    Returns
    -------
    dict
        mzQC-structured dictionary with quality metrics.
    """
    spectra = exp.getSpectra()
    ms1_count = 0
    ms2_count = 0
    tic_values = []
    rt_values = []
    precursor_mzs = []

    for spec in spectra:
        level = spec.getMSLevel()
        rt = spec.getRT()
        rt_values.append(rt)
        _, intensities = spec.get_peaks()
        tic = float(intensities.sum()) if len(intensities) > 0 else 0.0

        if level == 1:
            ms1_count += 1
            tic_values.append(tic)
        elif level == 2:
            ms2_count += 1
            for prec in spec.getPrecursors():
                precursor_mzs.append(prec.getMZ())

    tic_mean = sum(tic_values) / len(tic_values) if tic_values else 0.0
    tic_std = (
        math.sqrt(sum((t - tic_mean) ** 2 for t in tic_values) / len(tic_values))
        if tic_values
        else 0.0
    )

    rt_range = [min(rt_values), max(rt_values)] if rt_values else [0.0, 0.0]

    metrics = [
        {
            "accession": "QC:0000005",
            "name": "Number of MS1 spectra",
            "value": ms1_count,
        },
        {
            "accession": "QC:0000006",
            "name": "Number of MS2 spectra",
            "value": ms2_count,
        },
        {
            "accession": "QC:0000048",
            "name": "Mean TIC",
            "value": round(tic_mean, 4),
        },
        {
            "accession": "QC:0000049",
            "name": "TIC standard deviation",
            "value": round(tic_std, 4),
        },
        {
            "accession": "QC:0000019",
            "name": "RT range (seconds)",
            "value": [round(v, 2) for v in rt_range],
        },
        {
            "accession": "QC:0000029",
            "name": "Number of unique precursors",
            "value": len(set(round(mz, 4) for mz in precursor_mzs)),
        },
    ]

    mzqc = {
        "mzQC": {
            "version": "1.0.0",
            "creationDate": datetime.now(timezone.utc).isoformat(),
            "runQualities": [
                {
                    "metadata": {
                        "inputFiles": [{"name": input_file, "fileFormat": "mzML"}],
                    },
                    "qualityMetrics": metrics,
                }
            ],
        }
    }
    return mzqc


def main():
    parser = argparse.ArgumentParser(
        description="Generate mzQC JSON from an mzML file."
    )
    parser.add_argument("--input", required=True, metavar="FILE", help="Path to mzML file")
    parser.add_argument("--output", required=True, metavar="FILE", help="Output mzQC JSON path")
    args = parser.parse_args()

    exp = oms.MSExperiment()
    oms.MzMLFile().load(args.input, exp)

    mzqc = generate_mzqc(exp, input_file=args.input)

    with open(args.output, "w") as fh:
        json.dump(mzqc, fh, indent=2)

    n_metrics = len(mzqc["mzQC"]["runQualities"][0]["qualityMetrics"])
    print(f"mzQC written to {args.output} ({n_metrics} metrics)")


if __name__ == "__main__":
    main()
