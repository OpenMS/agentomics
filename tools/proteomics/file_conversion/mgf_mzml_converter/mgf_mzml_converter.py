"""
MGF ↔ mzML Converter
=====================
Bidirectional conversion between MGF (Mascot Generic Format) and mzML.
Combines and improves upon the separate mgf_to_mzml_converter and
mzml_to_mgf_converter tools.

Conversion direction is inferred automatically from the input file extension
when ``--direction`` is not supplied:
    - ``.mgf``  → mzML
    - ``.mzML`` / ``.mzml`` → MGF

Features
--------
- MGF → mzML: preserves title, precursor m/z, charge, retention time
- mzML → MGF: full filter suite (charge, RT range, m/z range, min peaks,
  min base-peak intensity, MS level)
- Auto-detection of conversion direction from file extension
- Conversion statistics printed to stdout

Usage
-----
    # Auto-detect direction from extensions
    python mgf_mzml_converter.py --input spectra.mgf --output spectra.mzML
    python mgf_mzml_converter.py --input run.mzML --output spectra.mgf

    # Explicit direction
    python mgf_mzml_converter.py --direction mgf2mzml --input spectra.mgf --output out.mzML
    python mgf_mzml_converter.py --direction mzml2mgf --input run.mzML --output out.mgf

    # mzML → MGF with filters
    python mgf_mzml_converter.py --input run.mzML --output out.mgf \\
        --charge 2 3 --rt-min 600 --rt-max 1800 --mz-min 400 --mz-max 1200 \\
        --min-peaks 5 --min-intensity 1000
"""

import os
import sys
from typing import List, Optional

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")

PROTON = 1.007276


# ---------------------------------------------------------------------------
# MGF → mzML
# ---------------------------------------------------------------------------


def parse_mgf(input_path: str) -> List[dict]:
    """Parse an MGF file and return a list of spectrum dicts.

    Each dict has keys: title, pepmass, charge, rt, peaks (list of
    (mz, intensity) tuples).

    Parameters
    ----------
    input_path:
        Path to the MGF file.

    Returns
    -------
    list of dict
    """
    spectra = []
    current: Optional[dict] = None

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
                    charge_str = line[7:].replace("+", "").replace("-", "").strip()
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
                            current["peaks"].append((float(parts[0]), float(parts[1])))
                        except ValueError:
                            pass

    return spectra


def convert_mgf_to_mzml(input_path: str, output_path: str) -> dict:
    """Convert an MGF file to mzML format.

    Parameters
    ----------
    input_path:
        Path to the input MGF file.
    output_path:
        Path to the output mzML file.

    Returns
    -------
    dict
        ``{"spectra_converted": int}``
    """
    mgf_spectra = parse_mgf(input_path)
    exp = oms.MSExperiment()

    for i, spec_data in enumerate(mgf_spectra):
        spectrum = oms.MSSpectrum()
        spectrum.setMSLevel(2)
        spectrum.setRT(spec_data["rt"])
        spectrum.setNativeID(spec_data["title"] if spec_data["title"] else f"index={i}")

        if spec_data["pepmass"] > 0:
            prec = oms.Precursor()
            prec.setMZ(spec_data["pepmass"])
            if spec_data["charge"] > 0:
                prec.setCharge(spec_data["charge"])
            spectrum.setPrecursors([prec])

        if spec_data["peaks"]:
            mzs = [p[0] for p in spec_data["peaks"]]
            intensities = [p[1] for p in spec_data["peaks"]]
            spectrum.set_peaks((mzs, intensities))

        exp.addSpectrum(spectrum)

    oms.MzMLFile().store(output_path, exp)
    return {"spectra_converted": len(mgf_spectra)}


# ---------------------------------------------------------------------------
# mzML → MGF
# ---------------------------------------------------------------------------


def load_mzml(input_path: str) -> oms.MSExperiment:
    """Load an mzML file into an MSExperiment.

    Parameters
    ----------
    input_path:
        Path to the mzML file.

    Returns
    -------
    oms.MSExperiment
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)
    return exp


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
    """Check whether a spectrum passes all active filters.

    Parameters
    ----------
    spectrum:
        The spectrum to check.
    min_peaks:
        Minimum number of fragment peaks required.
    charges:
        Allowed charge states (precursor).  ``None`` means any charge.
    rt_min:
        Minimum retention time in seconds.
    rt_max:
        Maximum retention time in seconds.
    mz_min:
        Minimum precursor m/z.
    mz_max:
        Maximum precursor m/z.
    min_intensity:
        Minimum base-peak intensity.

    Returns
    -------
    bool
    """
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
    """Convert a single spectrum to an MGF block string.

    Parameters
    ----------
    spectrum:
        The spectrum to convert.
    index:
        Fallback index used in the title when the native ID is empty.

    Returns
    -------
    str
        MGF block ending with a blank line.
    """
    lines = ["BEGIN IONS"]

    native_id = spectrum.getNativeID() if spectrum.getNativeID() else f"index={index}"
    lines.append(f"TITLE={native_id}")
    lines.append(f"RTINSECONDS={spectrum.getRT():.4f}")

    precursors = spectrum.getPrecursors()
    if precursors:
        prec = precursors[0]
        lines.append(f"PEPMASS={prec.getMZ():.6f}")
        if prec.getCharge() > 0:
            lines.append(f"CHARGE={prec.getCharge()}+")

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
    """Convert an mzML file to MGF format with optional filtering.

    Parameters
    ----------
    input_path:
        Path to the input mzML file.
    output_path:
        Path to the output MGF file.
    ms_level:
        MS level of spectra to export (default: 2).
    min_peaks:
        Minimum fragment peaks required to include a spectrum.
    charges:
        Allowed precursor charge states.
    rt_min:
        Minimum retention time in seconds.
    rt_max:
        Maximum retention time in seconds.
    mz_min:
        Minimum precursor m/z.
    mz_max:
        Maximum precursor m/z.
    min_intensity:
        Minimum base-peak intensity.

    Returns
    -------
    dict
        ``{"total_spectra", "ms_level_spectra", "converted", "filtered_out"}``
    """
    exp = load_mzml(input_path)
    target_spectra = [s for s in exp if s.getMSLevel() == ms_level]

    converted = 0
    filtered_out = 0
    with open(output_path, "w") as fh:
        for i, spectrum in enumerate(target_spectra):
            if not passes_filters(spectrum, min_peaks, charges, rt_min, rt_max, mz_min, mz_max, min_intensity):
                filtered_out += 1
                continue
            fh.write(spectrum_to_mgf_block(spectrum, i) + "\n")
            converted += 1

    return {
        "total_spectra": exp.size(),
        "ms_level_spectra": len(target_spectra),
        "converted": converted,
        "filtered_out": filtered_out,
    }


# ---------------------------------------------------------------------------
# Synthetic test data helper
# ---------------------------------------------------------------------------


def create_synthetic_mzml(output_path: str, n_spectra: int = 5) -> None:
    """Create a synthetic mzML file with MS2 spectra for testing.

    Parameters
    ----------
    output_path:
        Path for the output mzML file.
    n_spectra:
        Number of MS2 spectra to generate (default: 5).
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


# ---------------------------------------------------------------------------
# Direction auto-detection
# ---------------------------------------------------------------------------


def _detect_direction(input_path: str, output_path: str) -> str:
    """Infer conversion direction from file extensions.

    Parameters
    ----------
    input_path:
        Path to the input file.
    output_path:
        Path to the output file.

    Returns
    -------
    str
        ``"mgf2mzml"`` or ``"mzml2mgf"``.

    Raises
    ------
    click.UsageError
        If the direction cannot be determined.
    """
    in_ext = os.path.splitext(input_path)[1].lower()
    if in_ext == ".mgf":
        return "mgf2mzml"
    if in_ext in (".mzml", ".mzxml"):
        return "mzml2mgf"
    out_ext = os.path.splitext(output_path)[1].lower()
    if out_ext in (".mzml", ".mzxml"):
        return "mgf2mzml"
    if out_ext == ".mgf":
        return "mzml2mgf"
    raise click.UsageError(
        "Cannot determine conversion direction from file extensions. "
        "Use --direction mgf2mzml or --direction mzml2mgf explicitly."
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


@click.command(help="Bidirectional MGF ↔ mzML converter with optional spectrum filtering.")
@click.option("--input", "input_path", required=True, help="Input file (MGF or mzML)")
@click.option("--output", "output_path", required=True, help="Output file (mzML or MGF)")
@click.option(
    "--direction",
    type=click.Choice(["mgf2mzml", "mzml2mgf"]),
    default=None,
    help="Conversion direction (auto-detected from extension if omitted).",
)
@click.option("--ms-level", type=int, default=2, help="MS level to export (mzml2mgf, default: 2)")
@click.option("--min-peaks", type=int, default=1, help="Min fragment peaks to keep (mzml2mgf, default: 1)")
@click.option(
    "--charge", multiple=True, type=int,
    help="Keep only these precursor charge states (mzml2mgf, repeatable).",
)
@click.option("--rt-min", type=float, default=None, help="Min retention time in seconds (mzml2mgf)")
@click.option("--rt-max", type=float, default=None, help="Max retention time in seconds (mzml2mgf)")
@click.option("--mz-min", type=float, default=None, help="Min precursor m/z (mzml2mgf)")
@click.option("--mz-max", type=float, default=None, help="Max precursor m/z (mzml2mgf)")
@click.option(
    "--min-intensity", type=float, default=None,
    help="Min base-peak intensity to keep a spectrum (mzml2mgf).",
)
def main(
    input_path, output_path, direction, ms_level, min_peaks, charge,
    rt_min, rt_max, mz_min, mz_max, min_intensity,
) -> None:
    """CLI entry point."""
    if direction is None:
        direction = _detect_direction(input_path, output_path)

    if direction == "mgf2mzml":
        stats = convert_mgf_to_mzml(input_path, output_path)
        print(f"Converted {stats['spectra_converted']} spectra  {input_path} → {output_path}")
    else:
        charges = tuple(charge) if charge else None
        stats = convert_mzml_to_mgf(
            input_path, output_path, ms_level, min_peaks, charges,
            rt_min, rt_max, mz_min, mz_max, min_intensity,
        )
        print(
            f"Converted {stats['converted']} / {stats['ms_level_spectra']} "
            f"MS{ms_level} spectra  {input_path} → {output_path} "
            f"({stats['filtered_out']} filtered out)"
        )


if __name__ == "__main__":
    main()
