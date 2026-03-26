"""Tests for dia_scorer."""

import csv
import os

import pytest

pyopenms = pytest.importorskip("pyopenms")
import numpy as np  # noqa: E402


def _make_dia_mzml(path, product_mzs=None, n_scans=10):
    """Create a synthetic DIA mzML with known product m/z peaks."""
    if product_mzs is None:
        product_mzs = [600.0, 650.0, 700.0]

    exp = pyopenms.MSExperiment()
    for i in range(n_scans):
        spec = pyopenms.MSSpectrum()
        spec.setMSLevel(2)
        spec.setRT(float(i) * 10.0)
        center = n_scans // 2
        scale = max(0.1, 1.0 - abs(i - center) / center)
        mzs = np.array(product_mzs, dtype=np.float64)
        ints = np.array([5000.0 * scale] * len(product_mzs), dtype=np.float64)
        spec.set_peaks((mzs, ints))
        exp.addSpectrum(spec)
    pyopenms.MzMLFile().store(path, exp)


def _make_transitions_tsv(path, precursor_mz=500.0, product_mzs=None):
    """Create a transitions TSV file."""
    if product_mzs is None:
        product_mzs = [600.0, 650.0, 700.0]

    rows = []
    for i, pmz in enumerate(product_mzs):
        rows.append({
            "PrecursorMz": str(precursor_mz),
            "ProductMz": str(pmz),
            "LibraryIntensity": str(100.0 - i * 10),
            "PeptideSequence": "PEPTIDEK",
            "ProteinName": "PROT1",
            "transition_name": f"tr_{i}",
            "transition_group_id": "pep_0",
        })

    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class TestDIAScorer:
    def test_load_transitions(self):
        import tempfile

        from dia_scorer import load_transitions

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("PrecursorMz\tProductMz\tLibraryIntensity\tPeptideSequence\t"
                    "ProteinName\ttransition_name\ttransition_group_id\n")
            f.write("500.0\t600.0\t100.0\tPEPTIDEK\tPROT1\ttr_0\tpep_0\n")
            path = f.name

        try:
            trans = load_transitions(path)
            assert len(trans) == 1
            assert trans[0]["PrecursorMz"] == 500.0
            assert trans[0]["ProductMz"] == 600.0
        finally:
            os.unlink(path)

    def test_score_dia_produces_output(self, tmp_path):
        from dia_scorer import score_dia

        mzml_path = str(tmp_path / "dia.mzML")
        tsv_path = str(tmp_path / "transitions.tsv")
        out_path = str(tmp_path / "scores.tsv")

        _make_dia_mzml(mzml_path)
        _make_transitions_tsv(tsv_path)

        n = score_dia(mzml_path, tsv_path, out_path)
        assert n >= 1
        assert os.path.exists(out_path)

    def test_score_output_has_columns(self, tmp_path):
        from dia_scorer import score_dia

        mzml_path = str(tmp_path / "dia.mzML")
        tsv_path = str(tmp_path / "transitions.tsv")
        out_path = str(tmp_path / "scores.tsv")

        _make_dia_mzml(mzml_path)
        _make_transitions_tsv(tsv_path)

        score_dia(mzml_path, tsv_path, out_path)

        with open(out_path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = list(reader)

        assert len(rows) >= 1
        assert "transition_group_id" in rows[0]
        assert "dotprod_score" in rows[0]
        assert "manhattan_score" in rows[0]
        assert "best_rt" in rows[0]
        assert "n_transitions" in rows[0]

    def test_score_dia_returns_int(self, tmp_path):
        from dia_scorer import score_dia

        mzml_path = str(tmp_path / "dia.mzML")
        tsv_path = str(tmp_path / "transitions.tsv")
        out_path = str(tmp_path / "scores.tsv")

        _make_dia_mzml(mzml_path)
        _make_transitions_tsv(tsv_path)

        result = score_dia(mzml_path, tsv_path, out_path)
        assert isinstance(result, int)

    def test_group_transitions(self):
        from dia_scorer import _group_transitions

        transitions = [
            {"transition_group_id": "g1", "ProductMz": 600.0},
            {"transition_group_id": "g1", "ProductMz": 650.0},
            {"transition_group_id": "g2", "ProductMz": 700.0},
        ]
        groups = _group_transitions(transitions)
        assert len(groups) == 2
        assert len(groups["g1"]) == 2
        assert len(groups["g2"]) == 1
