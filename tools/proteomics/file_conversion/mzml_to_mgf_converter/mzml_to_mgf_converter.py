"""
mzML to MGF Converter
=====================
Convert MS2 spectra from mzML format to MGF (Mascot Generic Format)
with optional filtering by charge state, retention time, precursor m/z,
and minimum peak count.

Usage
-----
    python mzml_to_mgf_converter.py --input run.mzML --output spectra.mgf
    python mzml_to_mgf_converter.py --input run.mzML --output spectra.mgf --charge 2 3
    python mzml_to_mgf_converter.py --input run.mzML --output spectra.mgf --rt-min 600 --rt-max 1800
    python mzml_to_mgf_converter.py --input run.mzML --output spectra.mgf --mz-min 400 --mz-max 1200
"""

import sys
from typing import List, Optional

import click

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


def passes_filters(
    spectrum: oms.MSSpectrum,
    min_peaks: int = 1,
    charges: Optional[tuple] = None,
    rt_min: Optional[float] = None,
    rt_max: Optional[float] = None,
    mz_min: Optional[float] = None,
    mz_max: Optional[float] = None,
    min_intensity: Optional[float] = None,
) -> bool:
    """Check whether a spectrum passes all active filters."""
    mz_array, intensity_array = spectrum.get_peaks()
    if len(mz_array) < min_peaks:
        return False

    if rt_min is not None and spectrum.getRT() < rt_min:
        return False
    if rt_max is not None and spectrum.getRT() > rt_max:
        return False

    precursors = spectrum.getPrecursors()
    if precursors:
        prec = precursors[0]
        if charges and prec.getCharge() not in charges:
            return False
        if mz_min is not None and prec.getMZ() < mz_min:
            return False
        if mz_max is not None and prec.getMZ() > mz_max:
            return False
    elif charges or mz_min is not None or mz_max is not None:
        return False

    if min_intensity is not None:
        if len(intensity_array) == 0 or max(intensity_array) < min_intensity:
            return False

    return True


def spectrum_to_mgf_block(spectrum: oms.MSSpectrum, index: int) -> str:
    """Convert a single spectrum to an MGF block string."""
    lines = ["BEGIN IONS"]

    native_id = spectrum.getNativeID() if spectrum.getNativeID() else f"index={index}"
    lines.append(f"TITLE={native_id}")

    rt = spectrum.getRT()
    lines.append(f"RTINSECONDS={rt:.4f}")

    precursors = spectrum.getPrecursors()
    if precursors:
        prec = precursors[0]
        mz = prec.getMZ()
        charge = prec.getCharge()
        lines.append(f"PEPMASS={mz:.6f}")
        if charge > 0:
            lines.append(f"CHARGE={charge}+")

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
    charges: Optional[tuple] = None,
    rt_min: Optional[float] = None,
    rt_max: Optional[float] = None,
    mz_min: Optional[float] = None,
    mz_max: Optional[float] = None,
    min_intensity: Optional[float] = None,
) -> dict:
    """Convert mzML to MGF format with optional filtering.

    Returns statistics about the conversion.
    """
    exp = load_mzml(input_path)
    spectra = get_spectra_by_level(exp, ms_level)

    converted = 0
    filtered_out = 0
    with open(output_path, "w") as fh:
        for i, spectrum in enumerate(spectra):
            if not passes_filters(
                spectrum, min_peaks, charges, rt_min, rt_max, mz_min, mz_max,
                min_intensity,
            ):
                filtered_out += 1
                continue
            block = spectrum_to_mgf_block(spectrum, i)
            fh.write(block + "\n")
            converted += 1

    return {
        "total_spectra": exp.size(),
        "ms_level_spectra": len(spectra),
        "converted": converted,
        "filtered_out": filtered_out,
    }


def create_synthetic_mzml(output_path: str, n_spectra: int = 5) -> None:
    """Create a synthetic mzML file with MS2 spectra for testing.

    Generates spectra with:
    - RT: 10.0, 10.5, 11.0, ... (0.5s apart)
    - precursor m/z: 500, 550, 600, ...
    - charge: alternating 2, 3
    - base peak intensity: 10000, 9000, 8000, ...
    """
    exp = oms.MSExperiment()

    ms1 = oms.MSSpectrum()
    ms1.setMSLevel(1)
    ms1.setRT(10.0)
    ms1.set_peaks(([100.0, 200.0, 300.0], [1000.0, 2000.0, 1500.0]))
    exp.addSpectrum(ms1)

    for i in range(n_spectra):
        ms2 = oms.MSSpectrum()
        ms2.setMSLevel(2)
        ms2.setRT(10.0 + i * 0.5)
        ms2.setNativeID(f"scan={i + 1}")

        prec = oms.Precursor()
        prec.setMZ(500.0 + i * 50.0)
        prec.setCharge(2 if i % 2 == 0 else 3)
        ms2.setPrecursors([prec])

        mzs = [100.0 + j * 50.0 for j in range(10)]
        intensities = [1000.0 * (10 - j) for j in range(10)]
        ms2.set_peaks((mzs, intensities))
        exp.addSpectrum(ms2)

    oms.MzMLFile().store(output_path, exp)


@click.command(help="Convert MS2 spectra from mzML to MGF format with optional filtering.")
@click.option("--input", "input", required=True, help="Input mzML file")
@click.option("--ms-level", type=int, default=2, help="MS level to extract (default: 2)")
@click.option("--min-peaks", type=int, default=1, help="Minimum peaks per spectrum (default: 1)")
@click.option("--charge", multiple=True, type=int, help="Keep only these charge states (repeatable)")
@click.option("--rt-min", type=float, default=None, help="Minimum retention time in seconds")
@click.option("--rt-max", type=float, default=None, help="Maximum retention time in seconds")
@click.option("--mz-min", type=float, default=None, help="Minimum precursor m/z")
@click.option("--mz-max", type=float, default=None, help="Maximum precursor m/z")
@click.option(
    "--min-intensity", type=float, default=None,
    help="Minimum base peak intensity to keep a spectrum",
)
@click.option("--output", required=True, help="Output MGF file")
def main(input, ms_level, min_peaks, charge, rt_min, rt_max, mz_min, mz_max,
         min_intensity, output) -> None:
    charges = tuple(charge) if charge else None
    stats = convert_mzml_to_mgf(
        input, output, ms_level, min_peaks, charges, rt_min, rt_max,
        mz_min, mz_max, min_intensity,
    )
    print(
        f"Converted {stats['converted']} / {stats['ms_level_spectra']} "
        f"MS{ms_level} spectra to {output} "
        f"({stats['filtered_out']} filtered out)"
    )


if __name__ == "__main__":
    main()
