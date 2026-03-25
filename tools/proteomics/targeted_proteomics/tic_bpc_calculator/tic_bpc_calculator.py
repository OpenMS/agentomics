"""
TIC/BPC Calculator
==================
Compute Total Ion Chromatogram (TIC) and Base Peak Chromatogram (BPC) from mzML files.

Features:
- TIC: sum of all peak intensities per scan
- BPC: maximum peak intensity per scan
- Filter by MS level
- TSV output

Usage
-----
    python tic_bpc_calculator.py --input run.mzML --ms-level 1
    python tic_bpc_calculator.py --input run.mzML --ms-level 1 --output chromatograms.tsv
"""

import csv
import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def compute_tic_bpc(
    input_path: str,
    ms_level: int = 1,
) -> list[dict]:
    """Compute TIC and BPC chromatograms from mzML.

    Parameters
    ----------
    input_path : str
        Path to mzML file.
    ms_level : int
        MS level to compute chromatograms for (default 1).

    Returns
    -------
    list[dict]
        List of dicts with keys: scan_index, rt, tic, bpc, bpc_mz.
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(input_path, exp)

    results = []
    for i in range(exp.getNrSpectra()):
        spec = exp.getSpectrum(i)
        if spec.getMSLevel() != ms_level:
            continue

        rt = spec.getRT()
        mzs, intensities = spec.get_peaks()

        tic = 0.0
        bpc = 0.0
        bpc_mz = 0.0

        if len(intensities) > 0:
            tic = float(sum(intensities))
            max_idx = 0
            for j in range(len(intensities)):
                if intensities[j] > bpc:
                    bpc = float(intensities[j])
                    max_idx = j
            if len(mzs) > 0:
                bpc_mz = float(mzs[max_idx])

        results.append({
            "scan_index": i,
            "rt": round(rt, 4),
            "tic": round(tic, 2),
            "bpc": round(bpc, 2),
            "bpc_mz": round(bpc_mz, 6),
        })

    return results


def create_synthetic_mzml(output_path: str, n_scans: int = 10) -> None:
    """Create a synthetic mzML file for testing.

    Parameters
    ----------
    output_path : str
        Path to write the synthetic mzML file.
    n_scans : int
        Number of MS1 scans to generate.
    """
    exp = oms.MSExperiment()

    for i in range(n_scans):
        spec = oms.MSSpectrum()
        spec.setMSLevel(1)
        spec.setRT(float(i) * 10.0)
        mzs = [100.0 + j * 100 for j in range(5)]
        ints = [1000.0 * (j + 1) for j in range(5)]
        spec.set_peaks((mzs, ints))
        exp.addSpectrum(spec)

    oms.MzMLFile().store(output_path, exp)


def write_tsv(results: list[dict], output_path: str) -> None:
    """Write TIC/BPC results to TSV file.

    Parameters
    ----------
    results : list[dict]
        List of chromatogram data points.
    output_path : str
        Path to output TSV file.
    """
    fieldnames = ["scan_index", "rt", "tic", "bpc", "bpc_mz"]
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


@click.command(help="Compute TIC and BPC chromatograms from mzML.")
@click.option("--input", "input", required=True, help="Path to input mzML file")
@click.option("--ms-level", type=int, default=1, help="MS level (default: 1)")
@click.option("--output", default=None, help="Output TSV file path")
def main(input, ms_level, output):
    results = compute_tic_bpc(input, ms_level)

    if output:
        write_tsv(results, output)
        print(f"Wrote {len(results)} chromatogram data points to {output}")
    else:
        print("scan_index\trt\ttic\tbpc\tbpc_mz")
        for r in results:
            print(f"{r['scan_index']}\t{r['rt']}\t{r['tic']}\t{r['bpc']}\t{r['bpc_mz']}")


if __name__ == "__main__":
    main()
