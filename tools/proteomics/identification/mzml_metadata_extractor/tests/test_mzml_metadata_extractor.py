"""Tests for mzml_metadata_extractor."""

import json

import numpy as np
import pytest

pytest.importorskip("pyopenms")


class TestMzmlMetadataExtractor:
    def _make_experiment(self, n_spectra=5):
        """Create a synthetic MSExperiment."""
        import pyopenms as oms

        exp = oms.MSExperiment()
        for i in range(n_spectra):
            spec = oms.MSSpectrum()
            spec.setMSLevel(1 if i % 2 == 0 else 2)
            spec.setRT(60.0 * i)
            mzs = np.array([100.0 + j for j in range(10)], dtype=np.float64)
            ints = np.array([1000.0 * (j + 1) for j in range(10)], dtype=np.float64)
            spec.set_peaks([mzs, ints])
            exp.addSpectrum(spec)
        return exp

    def test_extract_metadata(self):
        from mzml_metadata_extractor import extract_metadata

        exp = self._make_experiment(n_spectra=4)
        meta = extract_metadata(exp)
        assert meta["n_spectra"] == 4
        assert "1" in meta["ms_levels"] or "2" in meta["ms_levels"]

    def test_rt_range(self):
        from mzml_metadata_extractor import extract_metadata

        exp = self._make_experiment(n_spectra=3)
        meta = extract_metadata(exp)
        assert len(meta["rt_range_sec"]) == 2
        assert meta["rt_range_sec"][0] == 0.0

    def test_mz_range(self):
        from mzml_metadata_extractor import extract_metadata

        exp = self._make_experiment(n_spectra=2)
        meta = extract_metadata(exp)
        assert len(meta["mz_range"]) == 2
        assert meta["mz_range"][0] == 100.0

    def test_empty_experiment(self):
        import pyopenms as oms
        from mzml_metadata_extractor import extract_metadata

        exp = oms.MSExperiment()
        meta = extract_metadata(exp)
        assert meta["n_spectra"] == 0
        assert meta["rt_range_sec"] == []

    def test_format_metadata(self):
        from mzml_metadata_extractor import extract_metadata, format_metadata

        exp = self._make_experiment(n_spectra=2)
        meta = extract_metadata(exp)
        text = format_metadata(meta)
        assert "Total spectra" in text

    def test_write_json(self, tmp_path):
        from mzml_metadata_extractor import extract_metadata, write_json

        exp = self._make_experiment(n_spectra=2)
        meta = extract_metadata(exp)
        out = str(tmp_path / "meta.json")
        write_json(meta, out)
        with open(out) as fh:
            data = json.load(fh)
        assert data["n_spectra"] == 2
