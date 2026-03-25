"""Tests for neutral_loss_scanner."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestNeutralLossScanner:
    def test_scan_with_known_losses(self):
        from neutral_loss_scanner import create_synthetic_mzml, scan_neutral_losses

        losses = [97.977, 162.053]
        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "test.mzML")
            create_synthetic_mzml(mzml_path, precursor_mz=500.0, losses=losses)
            results = scan_neutral_losses(mzml_path, losses, tolerance=0.05)
            assert len(results) >= 2
            found_losses = {r["neutral_loss"] for r in results}
            for loss in losses:
                assert round(loss, 6) in found_losses

    def test_no_match(self):
        from neutral_loss_scanner import create_synthetic_mzml, scan_neutral_losses

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "test.mzML")
            create_synthetic_mzml(mzml_path, precursor_mz=500.0, losses=[97.977])
            results = scan_neutral_losses(mzml_path, [999.0], tolerance=0.02)
            assert len(results) == 0

    def test_result_keys(self):
        from neutral_loss_scanner import create_synthetic_mzml, scan_neutral_losses

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "test.mzML")
            create_synthetic_mzml(mzml_path, precursor_mz=500.0, losses=[97.977])
            results = scan_neutral_losses(mzml_path, [97.977], tolerance=0.05)
            assert len(results) > 0
            for r in results:
                assert "scan_index" in r
                assert "precursor_mz" in r
                assert "neutral_loss" in r
                assert "fragment_mz" in r

    def test_write_tsv(self):
        from neutral_loss_scanner import write_tsv

        results = [{"scan_index": 0, "rt": 10.0, "precursor_mz": 500.0,
                     "neutral_loss": 97.977, "fragment_mz": 402.023, "intensity": 5000.0, "delta_da": 0.0}]
        with tempfile.TemporaryDirectory() as tmpdir:
            out = os.path.join(tmpdir, "matches.tsv")
            write_tsv(results, out)
            assert os.path.exists(out)

    def test_create_synthetic_mzml(self):
        import pyopenms as oms
        from neutral_loss_scanner import create_synthetic_mzml

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "synthetic.mzML")
            create_synthetic_mzml(mzml_path)
            exp = oms.MSExperiment()
            oms.MzMLFile().load(mzml_path, exp)
            assert exp.getNrSpectra() >= 2
