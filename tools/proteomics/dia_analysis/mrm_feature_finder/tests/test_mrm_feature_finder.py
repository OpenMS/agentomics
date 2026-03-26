"""Tests for mrm_feature_finder."""

import csv
import os

import pytest

pyopenms = pytest.importorskip("pyopenms")
import numpy as np  # noqa: E402


def _make_chromatograms_mzml(path, native_ids, peak_rts, peak_width=5.0):
    """Create synthetic chromatograms with Gaussian peaks."""
    exp = pyopenms.MSExperiment()
    times = np.linspace(0, 100, 100)

    for native_id, peak_rt in zip(native_ids, peak_rts):
        chrom = pyopenms.MSChromatogram()
        chrom.setNativeID(native_id.encode())
        ints = np.array(
            [1000.0 * np.exp(-0.5 * ((t - peak_rt) / peak_width) ** 2) for t in times]
        )
        chrom.set_peaks((times.tolist(), ints.tolist()))
        exp.addChromatogram(chrom)

    pyopenms.MzMLFile().store(path, exp)


def _make_transitions_tsv(path, native_ids, precursor_mz=500.0, product_mzs=None):
    """Create a transitions TSV file matching chromatogram native IDs."""
    if product_mzs is None:
        product_mzs = [600.0 + i * 50 for i in range(len(native_ids))]

    rows = []
    for i, (nid, pmz) in enumerate(zip(native_ids, product_mzs)):
        rows.append({
            "PrecursorMz": str(precursor_mz),
            "ProductMz": str(pmz),
            "LibraryIntensity": str(100.0 - i * 10),
            "PeptideSequence": "PEPTIDEK",
            "ProteinName": "PROT1",
            "transition_name": nid,
            "transition_group_id": "pep_0",
        })

    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class TestMRMFeatureFinder:
    def test_find_single_feature(self, tmp_path):
        from mrm_feature_finder import find_mrm_features

        native_ids = ["tr_0"]
        chrom_path = str(tmp_path / "chromatograms.mzML")
        tsv_path = str(tmp_path / "transitions.tsv")
        out_path = str(tmp_path / "features.featureXML")

        _make_chromatograms_mzml(chrom_path, native_ids, [50.0])
        _make_transitions_tsv(tsv_path, native_ids)

        n = find_mrm_features(chrom_path, tsv_path, out_path)
        assert n >= 1
        assert os.path.exists(out_path)

    def test_feature_rt_near_peak(self, tmp_path):
        from mrm_feature_finder import find_mrm_features

        native_ids = ["tr_0"]
        chrom_path = str(tmp_path / "chromatograms.mzML")
        tsv_path = str(tmp_path / "transitions.tsv")
        out_path = str(tmp_path / "features.featureXML")

        _make_chromatograms_mzml(chrom_path, native_ids, [50.0])
        _make_transitions_tsv(tsv_path, native_ids)

        find_mrm_features(chrom_path, tsv_path, out_path)

        fm = pyopenms.FeatureMap()
        pyopenms.FeatureXMLFile().load(out_path, fm)
        assert fm.size() >= 1
        # Feature RT should be near the Gaussian peak center (50.0)
        feature_rt = fm[0].getRT()
        assert abs(feature_rt - 50.0) < 15.0

    def test_multiple_transitions(self, tmp_path):
        from mrm_feature_finder import find_mrm_features

        native_ids = ["tr_0", "tr_1", "tr_2"]
        chrom_path = str(tmp_path / "chromatograms.mzML")
        tsv_path = str(tmp_path / "transitions.tsv")
        out_path = str(tmp_path / "features.featureXML")

        _make_chromatograms_mzml(chrom_path, native_ids, [50.0, 50.0, 50.0])
        _make_transitions_tsv(tsv_path, native_ids)

        n = find_mrm_features(chrom_path, tsv_path, out_path)
        assert n >= 1

    def test_returns_int(self, tmp_path):
        from mrm_feature_finder import find_mrm_features

        native_ids = ["tr_0"]
        chrom_path = str(tmp_path / "chromatograms.mzML")
        tsv_path = str(tmp_path / "transitions.tsv")
        out_path = str(tmp_path / "features.featureXML")

        _make_chromatograms_mzml(chrom_path, native_ids, [50.0])
        _make_transitions_tsv(tsv_path, native_ids)

        result = find_mrm_features(chrom_path, tsv_path, out_path)
        assert isinstance(result, int)

    def test_load_transitions(self, tmp_path):
        from mrm_feature_finder import load_transitions_to_targeted

        tsv_path = str(tmp_path / "transitions.tsv")
        _make_transitions_tsv(tsv_path, ["tr_0", "tr_1"])

        te = load_transitions_to_targeted(tsv_path)
        assert len(te.getTransitions()) == 2
        assert len(te.getPeptides()) == 1
        assert len(te.getProteins()) == 1
