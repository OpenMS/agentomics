"""
mzML Metadata Extractor
========================
Extract instrument metadata from mzML files and output as JSON.

Extracts: instrument model, source file, software, contact, MS levels,
number of spectra, RT range, and scan settings.

Usage
-----
    python mzml_metadata_extractor.py --input run.mzML --output metadata.json
"""

import json

import click
import pyopenms as oms


def extract_metadata(exp: oms.MSExperiment) -> dict:
    """Extract metadata from a loaded MSExperiment.

    Parameters
    ----------
    exp:
        Loaded MSExperiment object.

    Returns
    -------
    dict
        Metadata dictionary.
    """
    spectra = exp.getSpectra()
    n_spectra = len(spectra)

    # Collect MS levels
    ms_levels = {}
    rt_min = float("inf")
    rt_max = float("-inf")
    mz_min = float("inf")
    mz_max = float("-inf")

    for spec in spectra:
        level = spec.getMSLevel()
        ms_levels[level] = ms_levels.get(level, 0) + 1
        rt = spec.getRT()
        if rt < rt_min:
            rt_min = rt
        if rt > rt_max:
            rt_max = rt
        mzs, _ = spec.get_peaks()
        if len(mzs) > 0:
            local_min = float(mzs.min())
            local_max = float(mzs.max())
            if local_min < mz_min:
                mz_min = local_min
            if local_max > mz_max:
                mz_max = local_max

    # Extract instrument info from ExperimentalSettings
    instrument = exp.getInstrument()
    instrument_name = instrument.getName() if instrument else ""
    instrument_vendor = instrument.getVendor() if instrument else ""
    instrument_model = instrument.getModel() if instrument else ""

    # Source files
    source_files = []
    for sf in exp.getSourceFiles():
        source_files.append({
            "name": sf.getNameOfFile(),
            "path": sf.getPathToFile(),
        })

    # Software
    software_list = []
    seen_sw = set()
    for spec in exp:
        for dp in spec.getDataProcessing():
            sw = dp.getSoftware()
            sw_key = (sw.getName(), sw.getVersion())
            if sw_key not in seen_sw:
                seen_sw.add(sw_key)
                software_list.append({
                    "name": sw.getName(),
                    "version": sw.getVersion(),
                })

    metadata = {
        "n_spectra": n_spectra,
        "ms_levels": {str(k): v for k, v in sorted(ms_levels.items())},
        "rt_range_sec": [round(rt_min, 4), round(rt_max, 4)] if n_spectra > 0 else [],
        "mz_range": [round(mz_min, 6), round(mz_max, 6)] if mz_min != float("inf") else [],
        "instrument": {
            "name": instrument_name,
            "vendor": instrument_vendor,
            "model": instrument_model,
        },
        "source_files": source_files,
        "software": software_list,
    }

    return metadata


def write_json(metadata: dict, output_path: str) -> None:
    """Write metadata to a JSON file.

    Parameters
    ----------
    metadata:
        Metadata dictionary.
    output_path:
        Output file path.
    """
    with open(output_path, "w") as fh:
        json.dump(metadata, fh, indent=2)


def format_metadata(metadata: dict) -> str:
    """Format metadata for console output.

    Parameters
    ----------
    metadata:
        Metadata dictionary.

    Returns
    -------
    str
        Formatted string.
    """
    lines = []
    lines.append(f"Total spectra     : {metadata['n_spectra']}")
    for level, count in metadata["ms_levels"].items():
        lines.append(f"  MS{level} spectra    : {count}")
    if metadata["rt_range_sec"]:
        rt_min, rt_max = metadata["rt_range_sec"]
        lines.append(f"RT range          : {rt_min:.2f} - {rt_max:.2f} s")
    if metadata["mz_range"]:
        mz_min, mz_max = metadata["mz_range"]
        lines.append(f"m/z range         : {mz_min:.4f} - {mz_max:.4f}")
    inst = metadata["instrument"]
    if inst["name"] or inst["vendor"] or inst["model"]:
        lines.append(f"Instrument        : {inst['name']} {inst['vendor']} {inst['model']}".strip())
    for sf in metadata["source_files"]:
        lines.append(f"Source file       : {sf['name']}")
    for sw in metadata["software"]:
        lines.append(f"Software          : {sw['name']} {sw['version']}")
    return "\n".join(lines)


@click.command(help="Extract instrument metadata from mzML files.")
@click.option("--input", "input", required=True, help="Input mzML file")
@click.option("--output", default=None, help="Output JSON file path")
def main(input, output):
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input, exp)

    metadata = extract_metadata(exp)
    print(format_metadata(metadata))

    if output:
        write_json(metadata, output)
        print(f"\nMetadata written to {output}")


if __name__ == "__main__":
    main()
