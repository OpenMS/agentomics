"""Tests for chromatogram_extractor."""

import csv
import os

import pytest

pyopenms = pytest.importorskip("pyopenms")


def _make_dia_mzml(path, target_mz=500.0, product_mz=600.0, n_scans=10):
    """Create a synthetic DIA mzML with known peaks."""
    import numpy as np

    exp = pyopenms.MSExperiment()
    for i in range(n_scans):
        spec = pyopenms.MSSpectrum()
        spec.setMSLevel(2)
        spec.setRT(float(i) * 10.0)
        center = n_scans // 2
        intensity = max(100.0, 10000.0 * max(0, 1 - abs(i - center) / center))
        mzs = np.array([product_mz - 50, product_mz, product_mz + 50], dtype=np.float64)
        ints = np.array([200.0, intensity, 150.0], dtype=np.float64)
        spec.set_peaks((mzs, ints))
        exp.addSpectrum(spec)
    pyopenms.MzMLFile().store(path, exp)


def _make_transitions_tsv(path, precursor_mz=500.0, product_mz=600.0):
    """Create a transitions TSV file."""
    rows = [
        {
            "PrecursorMz": str(precursor_mz),
            "ProductMz": str(product_mz),
            "LibraryIntensity": "100.0",
            "PeptideSequence": "PEPTIDEK",
            "ProteinName": "PROT1",
            "transition_name": "trans_1",
            "transition_group_id": "pep_1",
        },
    ]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class TestChromatogramExtractor:
    def test_load_transitions_tsv(self, tmp_path):
        from chromatogram_extractor import load_transitions_tsv

        tsv_path = str(tmp_path / "transitions.tsv")
        _make_transitions_tsv(tsv_path)

        te = load_transitions_tsv(tsv_path)
        assert len(te.getTransitions()) == 1
        assert len(te.getPeptides()) == 1
        assert len(te.getProteins()) == 1

        t = te.getTransitions()[0]
        assert abs(t.getPrecursorMZ() - 500.0) < 0.01
        assert abs(t.getProductMZ() - 600.0) < 0.01

    def test_extract_chromatograms(self, tmp_path):
        from chromatogram_extractor import extract_chromatograms

        mzml_path = str(tmp_path / "dia.mzML")
        tsv_path = str(tmp_path / "transitions.tsv")
        out_path = str(tmp_path / "chromatograms.mzML")

        _make_dia_mzml(mzml_path)
        _make_transitions_tsv(tsv_path)

        n = extract_chromatograms(mzml_path, tsv_path, out_path, mz_tol=10.0)
        assert n >= 1
        assert os.path.exists(out_path)

        # Verify the output mzML contains chromatograms
        out_exp = pyopenms.MSExperiment()
        pyopenms.MzMLFile().load(out_path, out_exp)
        assert out_exp.getNrChromatograms() >= 1

    def test_chromatogram_has_intensities(self, tmp_path):
        from chromatogram_extractor import extract_chromatograms

        mzml_path = str(tmp_path / "dia.mzML")
        tsv_path = str(tmp_path / "transitions.tsv")
        out_path = str(tmp_path / "chromatograms.mzML")

        _make_dia_mzml(mzml_path)
        _make_transitions_tsv(tsv_path)

        extract_chromatograms(mzml_path, tsv_path, out_path, mz_tol=10.0)

        out_exp = pyopenms.MSExperiment()
        pyopenms.MzMLFile().load(out_path, out_exp)
        chrom = out_exp.getChromatogram(0)
        times, intensities = chrom.get_peaks()
        assert len(times) > 0
        assert max(intensities) > 0

    def test_multiple_transitions(self, tmp_path):
        from chromatogram_extractor import extract_chromatograms

        mzml_path = str(tmp_path / "dia.mzML")
        tsv_path = str(tmp_path / "transitions.tsv")
        out_path = str(tmp_path / "chromatograms.mzML")

        import numpy as np

        exp = pyopenms.MSExperiment()
        for i in range(10):
            spec = pyopenms.MSSpectrum()
            spec.setMSLevel(2)
            spec.setRT(float(i) * 10.0)
            mzs = np.array([600.0, 650.0, 700.0], dtype=np.float64)
            ints = np.array([1000.0, 2000.0, 500.0], dtype=np.float64)
            spec.set_peaks((mzs, ints))
            exp.addSpectrum(spec)
        pyopenms.MzMLFile().store(mzml_path, exp)

        rows = [
            {"PrecursorMz": "500.0", "ProductMz": "600.0", "LibraryIntensity": "100.0",
             "PeptideSequence": "PEPTIDEK", "ProteinName": "PROT1",
             "transition_name": "tr_0", "transition_group_id": "pep_0"},
            {"PrecursorMz": "500.0", "ProductMz": "650.0", "LibraryIntensity": "90.0",
             "PeptideSequence": "PEPTIDEK", "ProteinName": "PROT1",
             "transition_name": "tr_1", "transition_group_id": "pep_0"},
        ]
        with open(tsv_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=rows[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

        n = extract_chromatograms(mzml_path, tsv_path, out_path, mz_tol=10.0)
        assert n == 2

    def test_returns_int(self, tmp_path):
        from chromatogram_extractor import extract_chromatograms

        mzml_path = str(tmp_path / "dia.mzML")
        tsv_path = str(tmp_path / "transitions.tsv")
        out_path = str(tmp_path / "chromatograms.mzML")

        _make_dia_mzml(mzml_path)
        _make_transitions_tsv(tsv_path)

        result = extract_chromatograms(mzml_path, tsv_path, out_path)
        assert isinstance(result, int)
