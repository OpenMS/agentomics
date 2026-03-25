"""Tests for psm_feature_extractor."""

import os
import tempfile

from conftest import requires_pyopenms

PROTON = 1.007276


def _create_test_data(tmp_dir):
    """Create test mzML and PSM TSV matching specific peptides."""
    import pyopenms as oms

    sequences = ["ACDEFGHIK", "MNPQRSTWY"]
    exp = oms.MSExperiment()

    for i, seq in enumerate(sequences):
        aa_seq = oms.AASequence.fromString(seq)
        mass = aa_seq.getMonoWeight()
        precursor_mz = (mass + 2 * PROTON) / 2

        # Generate theoretical spectrum to use as "experimental"
        tsg = oms.TheoreticalSpectrumGenerator()
        params = tsg.getParameters()
        params.setValue("add_b_ions", "true")
        params.setValue("add_y_ions", "true")
        tsg.setParameters(params)

        theo = oms.MSSpectrum()
        tsg.getSpectrum(theo, aa_seq, 1, 2)

        # Use theoretical peaks as experimental (perfect match scenario)
        ms2 = oms.MSSpectrum()
        ms2.setMSLevel(2)
        ms2.setRT(100.0 + i * 20)
        prec = oms.Precursor()
        prec.setMZ(precursor_mz)
        prec.setCharge(2)
        ms2.setPrecursors([prec])

        mzs, ints = theo.get_peaks()
        ms2.set_peaks((list(mzs), list(ints)))
        exp.addSpectrum(ms2)

    mzml_path = os.path.join(tmp_dir, "test.mzML")
    oms.MzMLFile().store(mzml_path, exp)

    # Create PSMs TSV
    psms_path = os.path.join(tmp_dir, "psms.tsv")
    with open(psms_path, "w") as fh:
        fh.write("sequence\tcharge\trt\tscan_index\n")
        fh.write("ACDEFGHIK\t2\t100.0\t0\n")
        fh.write("MNPQRSTWY\t2\t120.0\t1\n")

    return mzml_path, psms_path


@requires_pyopenms
def test_generate_theoretical_spectrum():
    from psm_feature_extractor import generate_theoretical_spectrum

    theo = generate_theoretical_spectrum("ACDEFGHIK", 2)
    assert theo.size() > 0


@requires_pyopenms
def test_compute_features_perfect_match():
    from psm_feature_extractor import compute_features, generate_theoretical_spectrum

    theo = generate_theoretical_spectrum("ACDEFGHIK", 2)
    # Use same spectrum as experimental
    features = compute_features(theo, theo, "ACDEFGHIK", 2, tolerance=0.02)
    assert features["matched_fraction"] == 1.0
    assert features["matched_ions"] == features["total_theoretical"]
    assert features["sequence_length"] == 9


@requires_pyopenms
def test_extract_features():
    from psm_feature_extractor import extract_features

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path, psms_path = _create_test_data(tmp)
        output_path = os.path.join(tmp, "features.tsv")

        stats = extract_features(mzml_path, psms_path, output_path)
        assert stats["total_psms"] == 2
        assert stats["matched"] == 2

        with open(output_path) as fh:
            lines = fh.readlines()
        assert lines[0].strip().startswith("sequence")
        assert len(lines) == 3  # header + 2 PSMs


@requires_pyopenms
def test_extract_features_content():
    import csv

    from psm_feature_extractor import extract_features

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path, psms_path = _create_test_data(tmp)
        output_path = os.path.join(tmp, "features.tsv")

        extract_features(mzml_path, psms_path, output_path)

        with open(output_path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = list(reader)

        assert len(rows) == 2
        # Perfect match: matched_fraction should be high
        for row in rows:
            assert float(row["matched_fraction"]) > 0.5
            assert int(row["matched_ions"]) > 0
