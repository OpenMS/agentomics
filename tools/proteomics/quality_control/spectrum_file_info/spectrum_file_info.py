"""
Mass Spectrum File Info
=======================
Read an mzML (or mzXML) file and print a summary of its contents:
number of spectra, MS levels, retention time range, m/z range,
and basic TIC statistics.

Usage
-----
    python spectrum_file_info.py --input sample.mzML
    python spectrum_file_info.py --input sample.mzML --tic
"""


import click
import pyopenms as oms


def summarise_experiment(exp: oms.MSExperiment) -> dict:
    """Summarise a loaded MSExperiment object.

    Parameters
    ----------
    exp:
        Loaded ``pyopenms.MSExperiment`` instance.

    Returns
    -------
    dict
        Summary statistics for the experiment.
    """
    spectra = exp.getSpectra()
    if not spectra:
        return {"n_spectra": 0}

    ms_levels = {}
    rt_min = float("inf")
    rt_max = float("-inf")
    mz_min = float("inf")
    mz_max = float("-inf")
    tic_values = []

    for spec in spectra:
        level = spec.getMSLevel()
        ms_levels[level] = ms_levels.get(level, 0) + 1
        rt = spec.getRT()
        rt_min = min(rt_min, rt)
        rt_max = max(rt_max, rt)

        mzs, intensities = spec.get_peaks()
        if len(mzs) > 0:
            mz_min = min(mz_min, float(mzs.min()))
            mz_max = max(mz_max, float(mzs.max()))
            tic_values.append(float(intensities.sum()))

    return {
        "n_spectra": len(spectra),
        "ms_levels": ms_levels,
        "rt_range_sec": (rt_min, rt_max),
        "mz_range": (mz_min, mz_max),
        "tic_total": sum(tic_values),
        "tic_max": max(tic_values) if tic_values else 0.0,
        "tic_per_spectrum": tic_values,
    }


def load_file(path: str) -> oms.MSExperiment:
    """Load an mzML file into an MSExperiment.

    Parameters
    ----------
    path:
        Path to the mzML file.

    Returns
    -------
    pyopenms.MSExperiment
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(path, exp)
    return exp


@click.command(help="Summarise an mzML file using pyopenms.")
@click.option("--input", "input", required=True, help="Path to an mzML file")
@click.option("--tic", is_flag=True, help="Print per-spectrum TIC values")
def main(input, tic):
    print(f"Loading {input} ...")
    exp = load_file(input)
    summary = summarise_experiment(exp)

    if summary["n_spectra"] == 0:
        print("No spectra found in file.")
        return

    print(f"\n{'File':<22}: {input}")
    print(f"{'Total spectra':<22}: {summary['n_spectra']}")
    for level, count in sorted(summary["ms_levels"].items()):
        print(f"  {'MS' + str(level) + ' spectra':<20}: {count}")
    rt_min, rt_max = summary["rt_range_sec"]
    print(
        f"{'RT range':<22}: {rt_min:.2f} – {rt_max:.2f} s  "
        f"({rt_min/60:.2f} – {rt_max/60:.2f} min)"
    )
    mz_lo, mz_hi = summary["mz_range"]
    print(f"{'m/z range':<22}: {mz_lo:.4f} – {mz_hi:.4f}")
    print(f"{'Total TIC':<22}: {summary['tic_total']:.3e}")
    print(f"{'Max spectrum TIC':<22}: {summary['tic_max']:.3e}")

    if tic:
        print("\n--- Per-spectrum TIC ---")
        for i, tic in enumerate(summary["tic_per_spectrum"], 1):
            print(f"  Spectrum {i:>5}: {tic:.3e}")


if __name__ == "__main__":
    main()
