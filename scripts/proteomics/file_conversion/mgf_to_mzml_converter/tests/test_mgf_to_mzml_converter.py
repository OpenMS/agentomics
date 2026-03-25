"""Tests for mgf_to_mzml_converter."""

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


@requires_pyopenms
def test_parse_mgf():
    from mgf_to_mzml_converter import parse_mgf

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
    from mgf_to_mzml_converter import convert_mgf_to_mzml

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

        # Check first spectrum
        s = exp[0]
        assert s.getMSLevel() == 2
        precursors = s.getPrecursors()
        assert len(precursors) == 1
        assert abs(precursors[0].getMZ() - 500.25) < 0.01
        assert precursors[0].getCharge() == 2

        mz_arr, int_arr = s.get_peaks()
        assert len(mz_arr) == 3


@requires_pyopenms
def test_roundtrip_mgf_mzml_mgf():
    """Test MGF -> mzML -> MGF roundtrip."""
    import sys

    from mgf_to_mzml_converter import convert_mgf_to_mzml
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "mzml_to_mgf_converter"))

    with tempfile.TemporaryDirectory() as tmp:
        mgf_path = os.path.join(tmp, "test.mgf")
        mzml_path = os.path.join(tmp, "test.mzML")

        with open(mgf_path, "w") as fh:
            fh.write(SAMPLE_MGF)

        stats = convert_mgf_to_mzml(mgf_path, mzml_path)
        assert stats["spectra_converted"] == 2

        # Verify mzML is valid
        import pyopenms as oms
        exp = oms.MSExperiment()
        oms.MzMLFile().load(mzml_path, exp)
        assert exp.size() == 2
