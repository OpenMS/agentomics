"""Tests for spectral_library_builder."""

import os
import tempfile

from conftest import requires_pyopenms

PROTON = 1.007276


def _create_test_data(tmp_dir):
    """Create test mzML and peptides TSV."""
    import pyopenms as oms

    # Create mzML
    exp = oms.MSExperiment()
    sequences = ["ACDEFGHIK", "MNPQRSTWY"]

    for i, seq in enumerate(sequences):
        aa_seq = oms.AASequence.fromString(seq)
        mass = aa_seq.getMonoWeight()
        precursor_mz = (mass + 2 * PROTON) / 2

        ms2 = oms.MSSpectrum()
        ms2.setMSLevel(2)
        ms2.setRT(100.0 + i * 20)
        prec = oms.Precursor()
        prec.setMZ(precursor_mz)
        prec.setCharge(2)
        ms2.setPrecursors([prec])
        ms2.set_peaks(([100.0 + j * 50 for j in range(8)], [1000.0 - j * 100 for j in range(8)]))
        exp.addSpectrum(ms2)

    mzml_path = os.path.join(tmp_dir, "test.mzML")
    oms.MzMLFile().store(mzml_path, exp)

    # Create peptides TSV
    peptides_path = os.path.join(tmp_dir, "peptides.tsv")
    with open(peptides_path, "w") as fh:
        fh.write("sequence\tcharge\trt\n")
        fh.write("ACDEFGHIK\t2\t100.0\n")
        fh.write("MNPQRSTWY\t2\t120.0\n")

    return mzml_path, peptides_path


@requires_pyopenms
def test_load_peptides_tsv():
    from spectral_library_builder import load_peptides_tsv

    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "peptides.tsv")
        with open(path, "w") as fh:
            fh.write("sequence\tcharge\trt\n")
            fh.write("ACDEFGHIK\t2\t100.0\n")

        peptides = load_peptides_tsv(path)
        assert len(peptides) == 1
        assert peptides[0]["sequence"] == "ACDEFGHIK"
        assert peptides[0]["charge"] == 2
        # m/z should be auto-calculated
        assert peptides[0]["mz"] > 0


@requires_pyopenms
def test_build_library():
    from spectral_library_builder import build_library

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path, peptides_path = _create_test_data(tmp)
        output_path = os.path.join(tmp, "library.msp")

        stats = build_library(mzml_path, peptides_path, output_path)
        assert stats["total_peptides"] == 2
        assert stats["matched_spectra"] == 2

        with open(output_path) as fh:
            content = fh.read()
        assert "Name: ACDEFGHIK/2" in content
        assert "Num peaks:" in content


@requires_pyopenms
def test_no_match():
    import pyopenms as oms
    from spectral_library_builder import build_library

    with tempfile.TemporaryDirectory() as tmp:
        # Create mzML with no MS2 spectra
        exp = oms.MSExperiment()
        ms1 = oms.MSSpectrum()
        ms1.setMSLevel(1)
        ms1.set_peaks(([100.0], [1000.0]))
        exp.addSpectrum(ms1)
        mzml_path = os.path.join(tmp, "test.mzML")
        oms.MzMLFile().store(mzml_path, exp)

        peptides_path = os.path.join(tmp, "peptides.tsv")
        with open(peptides_path, "w") as fh:
            fh.write("sequence\tcharge\trt\n")
            fh.write("ACDEFGHIK\t2\t100.0\n")

        output_path = os.path.join(tmp, "library.msp")
        stats = build_library(mzml_path, peptides_path, output_path)
        assert stats["matched_spectra"] == 0
