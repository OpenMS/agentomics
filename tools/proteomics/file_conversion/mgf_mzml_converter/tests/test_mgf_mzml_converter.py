"""Tests for mgf_mzml_converter."""

import os
import tempfile

from conftest import requires_pyopenms

SAMPLE_MGF = """\
BEGIN IONS
TITLE=scan=1
RTINSECONDS=120.5
PEPMASS=500.250000
CHARGE=2+
100.100000 1000.0000
200.200000 2000.0000
300.300000 1500.0000
END IONS

BEGIN IONS
TITLE=scan=2
RTINSECONDS=130.0
PEPMASS=600.350000
CHARGE=3+
150.150000 800.0000
250.250000 1200.0000
END IONS
"""


# ---------------------------------------------------------------------------
# MGF → mzML tests
# ---------------------------------------------------------------------------


@requires_pyopenms
def test_parse_mgf():
    from mgf_mzml_converter import parse_mgf

    with tempfile.TemporaryDirectory() as tmp:
        mgf_path = os.path.join(tmp, "test.mgf")
        with open(mgf_path, "w") as fh:
            fh.write(SAMPLE_MGF)

        spectra = parse_mgf(mgf_path)
        assert len(spectra) == 2
        assert spectra[0]["title"] == "scan=1"
        assert abs(spectra[0]["pepmass"] - 500.25) < 0.01
        assert spectra[0]["charge"] == 2
        assert len(spectra[0]["peaks"]) == 3
        assert len(spectra[1]["peaks"]) == 2


@requires_pyopenms
def test_convert_mgf_to_mzml():
    import pyopenms as oms
    from mgf_mzml_converter import convert_mgf_to_mzml

    with tempfile.TemporaryDirectory() as tmp:
        mgf_path = os.path.join(tmp, "test.mgf")
        mzml_path = os.path.join(tmp, "test.mzML")

        with open(mgf_path, "w") as fh:
            fh.write(SAMPLE_MGF)

        stats = convert_mgf_to_mzml(mgf_path, mzml_path)
        assert stats["spectra_converted"] == 2

        exp = oms.MSExperiment()
        oms.MzMLFile().load(mzml_path, exp)
        assert exp.size() == 2

        s = exp[0]
        assert s.getMSLevel() == 2
        precursors = s.getPrecursors()
        assert len(precursors) == 1
        assert abs(precursors[0].getMZ() - 500.25) < 0.01
        assert precursors[0].getCharge() == 2

        mz_arr, _ = s.get_peaks()
        assert len(mz_arr) == 3


# ---------------------------------------------------------------------------
# mzML → MGF tests
# ---------------------------------------------------------------------------


@requires_pyopenms
def test_create_synthetic_mzml():
    import pyopenms as oms
    from mgf_mzml_converter import create_synthetic_mzml

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        create_synthetic_mzml(mzml_path, n_spectra=3)

        exp = oms.MSExperiment()
        oms.MzMLFile().load(mzml_path, exp)
        assert exp.size() == 4  # 1 MS1 + 3 MS2


@requires_pyopenms
def test_convert_mzml_to_mgf():
    from mgf_mzml_converter import convert_mzml_to_mgf, create_synthetic_mzml

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


@requires_pyopenms
def test_mgf_content_format():
    from mgf_mzml_converter import convert_mzml_to_mgf, create_synthetic_mzml

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


@requires_pyopenms
def test_min_peaks_filter():
    from mgf_mzml_converter import convert_mzml_to_mgf, create_synthetic_mzml

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        mgf_path = os.path.join(tmp, "test.mgf")

        create_synthetic_mzml(mzml_path, n_spectra=3)
        stats = convert_mzml_to_mgf(mzml_path, mgf_path, min_peaks=100)
        assert stats["converted"] == 0
        assert stats["filtered_out"] == 3


@requires_pyopenms
def test_charge_filter():
    from mgf_mzml_converter import convert_mzml_to_mgf, create_synthetic_mzml

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        mgf_path = os.path.join(tmp, "test.mgf")

        create_synthetic_mzml(mzml_path, n_spectra=6)
        stats = convert_mzml_to_mgf(mzml_path, mgf_path, charges=(2,))
        assert stats["converted"] == 3
        assert stats["filtered_out"] == 3

        with open(mgf_path) as fh:
            content = fh.read()
        assert "CHARGE=2+" in content
        assert "CHARGE=3+" not in content


@requires_pyopenms
def test_rt_range_filter():
    from mgf_mzml_converter import convert_mzml_to_mgf, create_synthetic_mzml

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        mgf_path = os.path.join(tmp, "test.mgf")

        create_synthetic_mzml(mzml_path, n_spectra=5)
        stats = convert_mzml_to_mgf(mzml_path, mgf_path, rt_min=10.5, rt_max=11.5)
        assert stats["converted"] == 3
        assert stats["filtered_out"] == 2


@requires_pyopenms
def test_mz_range_filter():
    from mgf_mzml_converter import convert_mzml_to_mgf, create_synthetic_mzml

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        mgf_path = os.path.join(tmp, "test.mgf")

        create_synthetic_mzml(mzml_path, n_spectra=5)
        stats = convert_mzml_to_mgf(mzml_path, mgf_path, mz_min=550.0, mz_max=650.0)
        assert stats["converted"] == 3
        assert stats["filtered_out"] == 2


@requires_pyopenms
def test_min_intensity_filter():
    from mgf_mzml_converter import convert_mzml_to_mgf, create_synthetic_mzml

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        mgf_path = os.path.join(tmp, "test.mgf")

        create_synthetic_mzml(mzml_path, n_spectra=3)
        stats = convert_mzml_to_mgf(mzml_path, mgf_path, min_intensity=10001.0)
        assert stats["converted"] == 0

        stats = convert_mzml_to_mgf(mzml_path, mgf_path, min_intensity=5000.0)
        assert stats["converted"] == 3


@requires_pyopenms
def test_combined_filters():
    from mgf_mzml_converter import convert_mzml_to_mgf, create_synthetic_mzml

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        mgf_path = os.path.join(tmp, "test.mgf")

        create_synthetic_mzml(mzml_path, n_spectra=6)
        stats = convert_mzml_to_mgf(
            mzml_path, mgf_path,
            charges=(2,), rt_min=10.0, rt_max=11.5, mz_min=500.0, mz_max=700.0,
        )
        assert stats["converted"] == 2


# ---------------------------------------------------------------------------
# Direction auto-detection tests
# ---------------------------------------------------------------------------


@requires_pyopenms
def test_auto_detect_mgf_to_mzml():
    from mgf_mzml_converter import _detect_direction

    direction = _detect_direction("spectra.mgf", "spectra.mzML")
    assert direction == "mgf2mzml"


@requires_pyopenms
def test_auto_detect_mzml_to_mgf():
    from mgf_mzml_converter import _detect_direction

    direction = _detect_direction("run.mzML", "out.mgf")
    assert direction == "mzml2mgf"


@requires_pyopenms
def test_auto_detect_mzml_lowercase():
    from mgf_mzml_converter import _detect_direction

    direction = _detect_direction("run.mzml", "out.mgf")
    assert direction == "mzml2mgf"


@requires_pyopenms
def test_auto_detect_from_output_extension():
    """When input has unknown extension, infer from output extension."""
    from mgf_mzml_converter import _detect_direction

    direction = _detect_direction("data.unknown", "out.mzML")
    assert direction == "mgf2mzml"


@requires_pyopenms
def test_roundtrip_mgf_mzml_mgf():
    """MGF → mzML → MGF round-trip preserves spectrum count."""
    import pyopenms as oms
    from mgf_mzml_converter import convert_mgf_to_mzml, convert_mzml_to_mgf

    with tempfile.TemporaryDirectory() as tmp:
        mgf_in = os.path.join(tmp, "in.mgf")
        mzml_path = os.path.join(tmp, "intermediate.mzML")
        mgf_out = os.path.join(tmp, "out.mgf")

        with open(mgf_in, "w") as fh:
            fh.write(SAMPLE_MGF)

        stats1 = convert_mgf_to_mzml(mgf_in, mzml_path)
        assert stats1["spectra_converted"] == 2

        exp = oms.MSExperiment()
        oms.MzMLFile().load(mzml_path, exp)
        assert exp.size() == 2

        stats2 = convert_mzml_to_mgf(mzml_path, mgf_out)
        assert stats2["converted"] == 2

        with open(mgf_out) as fh:
            content = fh.read()
        assert content.count("BEGIN IONS") == 2
