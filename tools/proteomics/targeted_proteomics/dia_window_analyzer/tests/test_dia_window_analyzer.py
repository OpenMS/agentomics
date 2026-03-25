"""Tests for dia_window_analyzer."""

import os
import tempfile

from conftest import requires_pyopenms


@requires_pyopenms
class TestDiaWindowAnalyzer:
    def test_analyze_windows(self):
        from dia_window_analyzer import analyze_dia_windows, create_synthetic_dia_mzml

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "dia.mzML")
            create_synthetic_dia_mzml(mzml_path, n_windows=5, window_width=25.0)
            results = analyze_dia_windows(mzml_path)
            assert len(results) == 5
            for r in results:
                assert r["window_width"] == 25.0
                assert r["scan_count"] == 3  # 3 cycles

    def test_result_keys(self):
        from dia_window_analyzer import analyze_dia_windows, create_synthetic_dia_mzml

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "dia.mzML")
            create_synthetic_dia_mzml(mzml_path, n_windows=2)
            results = analyze_dia_windows(mzml_path)
            for r in results:
                assert "window_center" in r
                assert "window_lower" in r
                assert "window_upper" in r
                assert "window_width" in r
                assert "scan_count" in r

    def test_window_boundaries(self):
        from dia_window_analyzer import analyze_dia_windows, create_synthetic_dia_mzml

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "dia.mzML")
            create_synthetic_dia_mzml(mzml_path, n_windows=3, window_width=20.0)
            results = analyze_dia_windows(mzml_path)
            for r in results:
                assert r["window_upper"] > r["window_lower"]
                assert abs(r["window_upper"] - r["window_lower"] - r["window_width"]) < 0.01

    def test_write_tsv(self):
        from dia_window_analyzer import analyze_dia_windows, create_synthetic_dia_mzml, write_tsv

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "dia.mzML")
            create_synthetic_dia_mzml(mzml_path)
            results = analyze_dia_windows(mzml_path)
            out = os.path.join(tmpdir, "windows.tsv")
            write_tsv(results, out)
            assert os.path.exists(out)
