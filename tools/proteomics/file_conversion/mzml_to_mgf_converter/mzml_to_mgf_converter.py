"""
mzML to MGF Converter
=====================
Convert MS2 spectra from mzML format to MGF (Mascot Generic Format).

Usage
-----
    python mzml_to_mgf_converter.py --input run.mzML --ms-level 2 --output spectra.mgf
"""

import argparse
import sys
from typing import List

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

PROTON = 1.007276


def load_mzml(input_path: str) -> oms.MSExperiment:
    """Load an mzML file into an MSExperiment."""
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)
    return exp


def get_spectra_by_level(exp: oms.MSExperiment, ms_level: int = 2) -> List[oms.MSSpectrum]:
    """Extract spectra of a given MS level from an experiment."""
    return [s for s in exp if s.getMSLevel() == ms_level]


def spectrum_to_mgf_block(spectrum: oms.MSSpectrum, index: int) -> str:
    """Convert a single spectrum to an MGF block string."""
    lines = ["BEGIN IONS"]

    # Title
    native_id = spectrum.getNativeID() if spectrum.getNativeID() else f"index={index}"
    lines.append(f"TITLE={native_id}")

    # Retention time
    rt = spectrum.getRT()
    lines.append(f"RTINSECONDS={rt:.4f}")

    # Precursor info
    precursors = spectrum.getPrecursors()
    if precursors:
        prec = precursors[0]
        mz = prec.getMZ()
        charge = prec.getCharge()
        lines.append(f"PEPMASS={mz:.6f}")
        if charge > 0:
            lines.append(f"CHARGE={charge}+")

    # Peaks
    mz_array, intensity_array = spectrum.get_peaks()
    for mz_val, intensity_val in zip(mz_array, intensity_array):
        lines.append(f"{mz_val:.6f} {intensity_val:.4f}")

    lines.append("END IONS")
    lines.append("")
    return "\n".join(lines)


def convert_mzml_to_mgf(
    input_path: str,
    output_path: str,
    ms_level: int = 2,
    min_peaks: int = 1,
) -> dict:
    """Convert mzML to MGF format.

    Returns statistics about the conversion.
    """
    exp = load_mzml(input_path)
    spectra = get_spectra_by_level(exp, ms_level)

    converted = 0
    with open(output_path, "w") as fh:
        for i, spectrum in enumerate(spectra):
            mz_array, _ = spectrum.get_peaks()
            if len(mz_array) < min_peaks:
                continue
            block = spectrum_to_mgf_block(spectrum, i)
            fh.write(block + "\n")
            converted += 1

    return {
        "total_spectra": exp.size(),
        "ms_level_spectra": len(spectra),
        "converted": converted,
    }


def create_synthetic_mzml(output_path: str, n_spectra: int = 5) -> None:
    """Create a synthetic mzML file with MS2 spectra for testing."""
    exp = oms.MSExperiment()

    # Add an MS1 spectrum
    ms1 = oms.MSSpectrum()
    ms1.setMSLevel(1)
    ms1.setRT(10.0)
    ms1.set_peaks(([100.0, 200.0, 300.0], [1000.0, 2000.0, 1500.0]))
    exp.addSpectrum(ms1)

    # Add MS2 spectra
    for i in range(n_spectra):
        ms2 = oms.MSSpectrum()
        ms2.setMSLevel(2)
        ms2.setRT(10.0 + i * 0.5)
        ms2.setNativeID(f"scan={i + 1}")

        prec = oms.Precursor()
        prec.setMZ(500.0 + i * 50.0)
        prec.setCharge(2)
        ms2.setPrecursors([prec])

        mzs = [100.0 + j * 50.0 for j in range(10)]
        intensities = [1000.0 * (10 - j) for j in range(10)]
        ms2.set_peaks((mzs, intensities))
        exp.addSpectrum(ms2)

    oms.MzMLFile().store(output_path, exp)


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert MS2 spectra from mzML to MGF format.")
    parser.add_argument("--input", required=True, help="Input mzML file")
    parser.add_argument("--ms-level", type=int, default=2, help="MS level to extract (default: 2)")
    parser.add_argument("--min-peaks", type=int, default=1, help="Minimum peaks per spectrum (default: 1)")
    parser.add_argument("--output", required=True, help="Output MGF file")
    args = parser.parse_args()

    stats = convert_mzml_to_mgf(args.input, args.output, args.ms_level, args.min_peaks)
    print(f"Converted {stats['converted']} / {stats['ms_level_spectra']} MS{args.ms_level} spectra to {args.output}")


if __name__ == "__main__":
    main()
