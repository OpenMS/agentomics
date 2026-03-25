"""
MGF to mzML Converter
=====================
Convert MGF (Mascot Generic Format) spectra to mzML format.

Usage
-----
    python mgf_to_mzml_converter.py --input spectra.mgf --output spectra.mzML
"""

import argparse
import sys
from typing import List

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def parse_mgf(input_path: str) -> List[dict]:
    """Parse an MGF file and return a list of spectrum dicts.

    Each dict has keys: title, pepmass, charge, rt, peaks (list of (mz, intensity) tuples).
    """
    spectra = []
    current = None

    with open(input_path) as fh:
        for line in fh:
            line = line.strip()
            if line == "BEGIN IONS":
                current = {"title": "", "pepmass": 0.0, "charge": 0, "rt": 0.0, "peaks": []}
            elif line == "END IONS":
                if current is not None:
                    spectra.append(current)
                current = None
            elif current is not None:
                if line.startswith("TITLE="):
                    current["title"] = line[6:]
                elif line.startswith("PEPMASS="):
                    parts = line[8:].split()
                    current["pepmass"] = float(parts[0])
                elif line.startswith("CHARGE="):
                    charge_str = line[7:].replace("+", "").replace("-", "")
                    try:
                        current["charge"] = int(charge_str)
                    except ValueError:
                        current["charge"] = 0
                elif line.startswith("RTINSECONDS="):
                    current["rt"] = float(line[12:])
                elif line and line[0].isdigit():
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            mz = float(parts[0])
                            intensity = float(parts[1])
                            current["peaks"].append((mz, intensity))
                        except ValueError:
                            pass

    return spectra


def convert_mgf_to_mzml(input_path: str, output_path: str) -> dict:
    """Convert MGF to mzML format.

    Returns statistics about the conversion.
    """
    mgf_spectra = parse_mgf(input_path)
    exp = oms.MSExperiment()

    for i, spec_data in enumerate(mgf_spectra):
        spectrum = oms.MSSpectrum()
        spectrum.setMSLevel(2)
        spectrum.setRT(spec_data["rt"])
        spectrum.setNativeID(spec_data["title"] if spec_data["title"] else f"index={i}")

        # Set precursor
        if spec_data["pepmass"] > 0:
            prec = oms.Precursor()
            prec.setMZ(spec_data["pepmass"])
            if spec_data["charge"] > 0:
                prec.setCharge(spec_data["charge"])
            spectrum.setPrecursors([prec])

        # Set peaks
        if spec_data["peaks"]:
            mzs = [p[0] for p in spec_data["peaks"]]
            intensities = [p[1] for p in spec_data["peaks"]]
            spectrum.set_peaks((mzs, intensities))

        exp.addSpectrum(spectrum)

    oms.MzMLFile().store(output_path, exp)

    return {"spectra_converted": len(mgf_spectra)}


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert MGF to mzML format.")
    parser.add_argument("--input", required=True, help="Input MGF file")
    parser.add_argument("--output", required=True, help="Output mzML file")
    args = parser.parse_args()

    stats = convert_mgf_to_mzml(args.input, args.output)
    print(f"Converted {stats['spectra_converted']} spectra to {args.output}")


if __name__ == "__main__":
    main()
