"""Tests for tic_bpc_calculator."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestTicBpcCalculator:
    def test_compute_tic_bpc(self):
        from tic_bpc_calculator import compute_tic_bpc, create_synthetic_mzml

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "test.mzML")
            create_synthetic_mzml(mzml_path, n_scans=5)
            results = compute_tic_bpc(mzml_path, ms_level=1)
            assert len(results) == 5
            for r in results:
                assert r["tic"] > 0
                assert r["bpc"] > 0

    def test_tic_equals_sum(self):
        from tic_bpc_calculator import compute_tic_bpc, create_synthetic_mzml

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "test.mzML")
            create_synthetic_mzml(mzml_path, n_scans=1)
            results = compute_tic_bpc(mzml_path, ms_level=1)
            assert len(results) == 1
            # TIC should be sum of [1000, 2000, 3000, 4000, 5000] = 15000
            assert results[0]["tic"] == 15000.0

    def test_bpc_is_max(self):
        from tic_bpc_calculator import compute_tic_bpc, create_synthetic_mzml

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "test.mzML")
            create_synthetic_mzml(mzml_path, n_scans=1)
            results = compute_tic_bpc(mzml_path, ms_level=1)
            assert results[0]["bpc"] == 5000.0

    def test_result_keys(self):
        from tic_bpc_calculator import compute_tic_bpc, create_synthetic_mzml

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "test.mzML")
            create_synthetic_mzml(mzml_path, n_scans=1)
            results = compute_tic_bpc(mzml_path)
            for r in results:
                assert "scan_index" in r
                assert "rt" in r
                assert "tic" in r
                assert "bpc" in r
                assert "bpc_mz" in r

    def test_write_tsv(self):
        from tic_bpc_calculator import compute_tic_bpc, create_synthetic_mzml, write_tsv

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "test.mzML")
            create_synthetic_mzml(mzml_path)
            results = compute_tic_bpc(mzml_path)
            out = os.path.join(tmpdir, "chrom.tsv")
            write_tsv(results, out)
            assert os.path.exists(out)
