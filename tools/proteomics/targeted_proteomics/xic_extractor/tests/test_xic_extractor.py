"""Tests for xic_extractor."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestXicExtractor:
    def test_extract_xic(self):
        from xic_extractor import create_synthetic_mzml, extract_xic

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "test.mzML")
            create_synthetic_mzml(mzml_path, target_mz=524.265, n_scans=10)
            results = extract_xic(mzml_path, 524.265, ppm=10.0)
            assert len(results) == 10
            # Middle scan should have highest intensity
            intensities = [r["intensity"] for r in results]
            assert max(intensities) > 0

    def test_xic_no_peaks(self):
        from xic_extractor import create_synthetic_mzml, extract_xic

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "test.mzML")
            create_synthetic_mzml(mzml_path, target_mz=524.265, n_scans=5)
            results = extract_xic(mzml_path, 999.999, ppm=1.0)
            assert len(results) == 5
            assert all(r["intensity"] == 0.0 for r in results)

    def test_result_keys(self):
        from xic_extractor import create_synthetic_mzml, extract_xic

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "test.mzML")
            create_synthetic_mzml(mzml_path)
            results = extract_xic(mzml_path, 524.265, ppm=10.0)
            for r in results:
                assert "rt" in r
                assert "intensity" in r
                assert "mz" in r

    def test_write_tsv(self):
        from xic_extractor import create_synthetic_mzml, extract_xic, write_tsv

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "test.mzML")
            create_synthetic_mzml(mzml_path)
            results = extract_xic(mzml_path, 524.265, ppm=10.0)
            out = os.path.join(tmpdir, "xic.tsv")
            write_tsv(results, out)
            assert os.path.exists(out)
            with open(out) as f:
                lines = f.readlines()
            assert len(lines) > 1

    def test_create_synthetic(self):
        import pyopenms as oms
        from xic_extractor import create_synthetic_mzml

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "synthetic.mzML")
            create_synthetic_mzml(mzml_path, n_scans=5)
            exp = oms.MSExperiment()
            oms.MzMLFile().load(mzml_path, exp)
            assert exp.getNrSpectra() == 5
