"""
DIA Window Analyzer
===================
Report DIA isolation window scheme from mzML metadata.
Extracts precursor isolation window information from MS2 spectra.

Features:
- Extract isolation window center and width
- Detect DIA window scheme
- Report window overlap
- TSV output

Usage
-----
    python dia_window_analyzer.py --input dia.mzML
    python dia_window_analyzer.py --input dia.mzML --output windows.tsv
"""

import csv

import click
import pyopenms as oms


def analyze_dia_windows(input_path: str) -> list[dict]:
    """Analyze DIA isolation windows from mzML file.

    Parameters
    ----------
    input_path : str
        Path to mzML file.

    Returns
    -------
    list[dict]
        List of unique DIA windows with keys: window_center, window_lower,
        window_upper, window_width, scan_count.
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    window_counts = {}

    for i in range(exp.getNrSpectra()):
        spec = exp.getSpectrum(i)
        if spec.getMSLevel() != 2:
            continue

        precursors = spec.getPrecursors()
        if not precursors:
            continue

        prec = precursors[0]
        center = prec.getMZ()
        lower = prec.getIsolationWindowLowerOffset()
        upper = prec.getIsolationWindowUpperOffset()

        key = (round(center, 4), round(lower, 4), round(upper, 4))
        window_counts[key] = window_counts.get(key, 0) + 1

    results = []
    for (center, lower, upper), count in sorted(window_counts.items()):
        results.append({
            "window_center": center,
            "window_lower": round(center - lower, 4),
            "window_upper": round(center + upper, 4),
            "window_width": round(lower + upper, 4),
            "scan_count": count,
        })

    return results


def create_synthetic_dia_mzml(output_path: str, n_windows: int = 5, window_width: float = 25.0) -> None:
    """Create a synthetic DIA mzML file for testing.

    Parameters
    ----------
    output_path : str
        Path to write the synthetic mzML file.
    n_windows : int
        Number of DIA windows.
    window_width : float
        Width of each window in Da.
    """
    exp = oms.MSExperiment()

    start_mz = 400.0
    for cycle in range(3):
        # MS1 survey scan
        ms1 = oms.MSSpectrum()
        ms1.setMSLevel(1)
        ms1.setRT(float(cycle) * 10.0)
        ms1.set_peaks(([500.0], [10000.0]))
        exp.addSpectrum(ms1)

        # DIA MS2 scans
        for w in range(n_windows):
            ms2 = oms.MSSpectrum()
            ms2.setMSLevel(2)
            ms2.setRT(float(cycle) * 10.0 + float(w) + 1.0)

            center = start_mz + w * window_width + window_width / 2
            prec = oms.Precursor()
            prec.setMZ(center)
            prec.setIsolationWindowLowerOffset(window_width / 2)
            prec.setIsolationWindowUpperOffset(window_width / 2)
            ms2.setPrecursors([prec])

            ms2.set_peaks(([center - 5, center, center + 5], [500.0, 1000.0, 300.0]))
            exp.addSpectrum(ms2)

    oms.MzMLFile().store(output_path, exp)


def write_tsv(results: list[dict], output_path: str) -> None:
    """Write DIA window results to TSV file.

    Parameters
    ----------
    results : list[dict]
        List of window result dictionaries.
    output_path : str
        Path to output TSV file.
    """
    fieldnames = ["window_center", "window_lower", "window_upper", "window_width", "scan_count"]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


@click.command(help="Report DIA isolation window scheme from mzML metadata.")
@click.option("--input", "input", required=True, help="Path to input mzML file")
@click.option("--output", default=None, help="Output TSV file path")
def main(input, output):
    results = analyze_dia_windows(input)

    if output:
        write_tsv(results, output)
        print(f"Wrote {len(results)} DIA windows to {output}")
    else:
        print("window_center\twindow_lower\twindow_upper\twindow_width\tscan_count")
        for r in results:
            print(
                f"{r['window_center']}\t{r['window_lower']}\t{r['window_upper']}\t"
                f"{r['window_width']}\t{r['scan_count']}"
            )


if __name__ == "__main__":
    main()
