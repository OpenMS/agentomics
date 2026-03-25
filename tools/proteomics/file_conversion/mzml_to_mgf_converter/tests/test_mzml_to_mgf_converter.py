"""Tests for mzml_to_mgf_converter."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


def test_create_synthetic_mzml():
    import pyopenms as oms
    from mzml_to_mgf_converter import create_synthetic_mzml

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        create_synthetic_mzml(mzml_path, n_spectra=3)

        exp = oms.MSExperiment()
        oms.MzMLFile().load(mzml_path, exp)
        assert exp.size() == 4  # 1 MS1 + 3 MS2


def test_convert_mzml_to_mgf():
    from mzml_to_mgf_converter import convert_mzml_to_mgf, create_synthetic_mzml

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        mgf_path = os.path.join(tmp, "test.mgf")

        create_synthetic_mzml(mzml_path, n_spectra=3)
        stats = convert_mzml_to_mgf(mzml_path, mgf_path)

        assert stats["converted"] == 3
        assert stats["ms_level_spectra"] == 3
        assert stats["filtered_out"] == 0

        with open(mgf_path) as fh:
            content = fh.read()
        assert content.count("BEGIN IONS") == 3
        assert "PEPMASS=" in content
        assert "CHARGE=" in content


def test_mgf_content_format():
    from mzml_to_mgf_converter import convert_mzml_to_mgf, create_synthetic_mzml

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        mgf_path = os.path.join(tmp, "test.mgf")

        create_synthetic_mzml(mzml_path, n_spectra=1)
        convert_mzml_to_mgf(mzml_path, mgf_path)

        with open(mgf_path) as fh:
            lines = fh.readlines()

        assert lines[0].strip() == "BEGIN IONS"
        assert any(line.startswith("TITLE=") for line in lines)
        assert any(line.startswith("PEPMASS=") for line in lines)


def test_min_peaks_filter():
    from mzml_to_mgf_converter import convert_mzml_to_mgf, create_synthetic_mzml

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        mgf_path = os.path.join(tmp, "test.mgf")

        create_synthetic_mzml(mzml_path, n_spectra=3)
        stats = convert_mzml_to_mgf(mzml_path, mgf_path, min_peaks=100)
        assert stats["converted"] == 0
        assert stats["filtered_out"] == 3


def test_charge_filter():
    """Synthetic spectra alternate charge 2, 3. Filter for charge 2 only."""
    from mzml_to_mgf_converter import convert_mzml_to_mgf, create_synthetic_mzml

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        mgf_path = os.path.join(tmp, "test.mgf")

        create_synthetic_mzml(mzml_path, n_spectra=6)
        stats = convert_mzml_to_mgf(mzml_path, mgf_path, charges=(2,))
        # Spectra 0, 2, 4 have charge 2; spectra 1, 3, 5 have charge 3
        assert stats["converted"] == 3
        assert stats["filtered_out"] == 3

        with open(mgf_path) as fh:
            content = fh.read()
        assert "CHARGE=2+" in content
        assert "CHARGE=3+" not in content


def test_charge_filter_multiple():
    """Filter for both charge 2 and 3 — should keep all."""
    from mzml_to_mgf_converter import convert_mzml_to_mgf, create_synthetic_mzml

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        mgf_path = os.path.join(tmp, "test.mgf")

        create_synthetic_mzml(mzml_path, n_spectra=4)
        stats = convert_mzml_to_mgf(mzml_path, mgf_path, charges=(2, 3))
        assert stats["converted"] == 4
        assert stats["filtered_out"] == 0


def test_rt_range_filter():
    """Synthetic spectra at RT 10.0, 10.5, 11.0, 11.5, 12.0. Filter RT 10.5-11.5."""
    from mzml_to_mgf_converter import convert_mzml_to_mgf, create_synthetic_mzml

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        mgf_path = os.path.join(tmp, "test.mgf")

        create_synthetic_mzml(mzml_path, n_spectra=5)
        stats = convert_mzml_to_mgf(
            mzml_path, mgf_path, rt_min=10.5, rt_max=11.5,
        )
        # RT 10.5, 11.0, 11.5 pass; RT 10.0, 12.0 filtered out
        assert stats["converted"] == 3
        assert stats["filtered_out"] == 2


def test_mz_range_filter():
    """Synthetic precursor m/z: 500, 550, 600, 650, 700. Filter 550-650."""
    from mzml_to_mgf_converter import convert_mzml_to_mgf, create_synthetic_mzml

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        mgf_path = os.path.join(tmp, "test.mgf")

        create_synthetic_mzml(mzml_path, n_spectra=5)
        stats = convert_mzml_to_mgf(
            mzml_path, mgf_path, mz_min=550.0, mz_max=650.0,
        )
        # m/z 550, 600, 650 pass; 500, 700 filtered out
        assert stats["converted"] == 3
        assert stats["filtered_out"] == 2


def test_min_intensity_filter():
    """All synthetic spectra have base peak intensity 10000. Filter > 10000 removes all."""
    from mzml_to_mgf_converter import convert_mzml_to_mgf, create_synthetic_mzml

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        mgf_path = os.path.join(tmp, "test.mgf")

        create_synthetic_mzml(mzml_path, n_spectra=3)
        stats = convert_mzml_to_mgf(mzml_path, mgf_path, min_intensity=10001.0)
        assert stats["converted"] == 0

        stats = convert_mzml_to_mgf(mzml_path, mgf_path, min_intensity=5000.0)
        assert stats["converted"] == 3


def test_combined_filters():
    """Apply charge + RT + m/z filters together."""
    from mzml_to_mgf_converter import convert_mzml_to_mgf, create_synthetic_mzml

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        mgf_path = os.path.join(tmp, "test.mgf")

        # 6 spectra: RT 10.0-12.5, m/z 500-750, charge alternating 2/3
        create_synthetic_mzml(mzml_path, n_spectra=6)
        stats = convert_mzml_to_mgf(
            mzml_path, mgf_path,
            charges=(2,), rt_min=10.0, rt_max=11.5, mz_min=500.0, mz_max=700.0,
        )
        # charge 2: indices 0, 2, 4 (RT 10.0/11.0/12.0, m/z 500/600/700)
        # RT filter: removes index 4 (RT 12.0)
        # m/z filter: all 500, 600, 700 pass (700 is at boundary, passes <=)
        # Result: indices 0, 2 pass
        assert stats["converted"] == 2
