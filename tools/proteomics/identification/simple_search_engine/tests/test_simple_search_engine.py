"""Tests for simple_search_engine."""

import pytest

pyopenms = pytest.importorskip("pyopenms")


def _make_synthetic_mzml(path, sequence="PEPTIDEK", charge=2):
    """Create a synthetic mzML with a theoretical MS2 spectrum for the given peptide."""
    import numpy as np
    import pyopenms as oms

    aa_seq = oms.AASequence.fromString(sequence)

    # Generate theoretical MS2 spectrum (b/y ions) with charge 1 fragments
    spec = oms.MSSpectrum()
    generator = oms.TheoreticalSpectrumGenerator()
    params = generator.getParameters()
    params.setValue("add_b_ions", "true")
    params.setValue("add_y_ions", "true")
    params.setValue("add_metainfo", "true")
    generator.setParameters(params)
    generator.getSpectrum(spec, aa_seq, 1, 1)
    spec.sortByPosition()

    # Set realistic intensities
    mzs, ints = spec.get_peaks()
    ints = np.ones_like(ints) * 10000.0
    spec.set_peaks([mzs, ints])

    # Set MS2 level, retention time, and native ID
    spec.setMSLevel(2)
    spec.setRT(100.0)
    spec.setNativeID("controllerType=0 controllerNumber=1 scan=2")

    precursor_mz = aa_seq.getMZ(charge)
    precursor = oms.Precursor()
    precursor.setMZ(precursor_mz)
    precursor.setCharge(charge)
    spec.setPrecursors([precursor])

    # Create an MS1 spectrum
    ms1 = oms.MSSpectrum()
    ms1.setMSLevel(1)
    ms1.setRT(99.0)
    ms1.setNativeID("controllerType=0 controllerNumber=1 scan=1")
    ms1.set_peaks([
        np.array(
            [precursor_mz - 1, precursor_mz, precursor_mz + 0.5],
            dtype=np.float64,
        ),
        np.array([5000.0, 100000.0, 50000.0], dtype=np.float64),
    ])

    exp = oms.MSExperiment()
    exp.addSpectrum(ms1)
    exp.addSpectrum(spec)
    oms.MzMLFile().store(str(path), exp)


def _make_fasta(path, proteins):
    """Create a FASTA file with the given proteins."""
    import pyopenms as oms

    entries = []
    for ident, seq in proteins:
        entry = oms.FASTAEntry()
        entry.identifier = ident
        entry.description = f"Test protein {ident}"
        entry.sequence = seq
        entries.append(entry)
    oms.FASTAFile().store(str(path), entries)


class TestSimpleSearchEngine:
    def test_search_finds_peptidek(self, tmp_path):
        from simple_search_engine import search

        mzml_path = tmp_path / "run.mzML"
        fasta_path = tmp_path / "proteins.fasta"
        output_path = tmp_path / "results.idXML"

        _make_synthetic_mzml(mzml_path, "PEPTIDEK", charge=2)
        # PEPTIDEK must be a proper tryptic peptide (at protein N-terminus)
        _make_fasta(fasta_path, [
            ("sp|P00001|PROT1", "PEPTIDEKAAAAR"),
        ])

        psm_count = search(
            str(mzml_path),
            str(fasta_path),
            str(output_path),
            precursor_tol=10.0,
            fragment_tol=10.0,
        )

        assert psm_count >= 1
        assert output_path.exists()

    def test_search_output_contains_peptidek(self, tmp_path):
        from simple_search_engine import search

        mzml_path = tmp_path / "run.mzML"
        fasta_path = tmp_path / "proteins.fasta"
        output_path = tmp_path / "results.idXML"

        _make_synthetic_mzml(mzml_path, "PEPTIDEK", charge=2)
        _make_fasta(fasta_path, [
            ("sp|P00001|PROT1", "PEPTIDEKAAAAR"),
        ])

        search(
            str(mzml_path),
            str(fasta_path),
            str(output_path),
            precursor_tol=10.0,
            fragment_tol=10.0,
        )

        # Load and verify idXML results
        import pyopenms as oms

        prot_ids = []
        pep_ids = oms.PeptideIdentificationList()
        oms.IdXMLFile().load(str(output_path), prot_ids, pep_ids)

        # Check that we got at least one peptide identification
        assert len(pep_ids) >= 1

        # Check that PEPTIDEK is among the identified sequences
        found_peptidek = False
        for pid in pep_ids:
            for hit in pid.getHits():
                seq_str = hit.getSequence().toString()
                if seq_str == "PEPTIDEK":
                    found_peptidek = True
                    break
        assert found_peptidek, "PEPTIDEK not found in search results"

    def test_search_returns_zero_for_no_match(self, tmp_path):
        from simple_search_engine import search

        mzml_path = tmp_path / "run.mzML"
        fasta_path = tmp_path / "proteins.fasta"
        output_path = tmp_path / "results.idXML"

        # Create mzML with PEPTIDEK spectrum but FASTA with unrelated protein
        _make_synthetic_mzml(mzml_path, "PEPTIDEK", charge=2)
        _make_fasta(fasta_path, [
            ("sp|P99999|UNRELATED", "MAAAAAAAAAAAAAAAAAAAAAAAAAAAAR"),
        ])

        psm_count = search(
            str(mzml_path),
            str(fasta_path),
            str(output_path),
            precursor_tol=10.0,
            fragment_tol=10.0,
        )

        # Should find no PSMs since the database has no matching peptide
        assert psm_count == 0

    def test_search_writes_valid_idxml(self, tmp_path):
        from simple_search_engine import search

        mzml_path = tmp_path / "run.mzML"
        fasta_path = tmp_path / "proteins.fasta"
        output_path = tmp_path / "results.idXML"

        _make_synthetic_mzml(mzml_path, "PEPTIDEK", charge=2)
        _make_fasta(fasta_path, [
            ("sp|P00001|PROT1", "PEPTIDEKAAAAR"),
        ])

        search(
            str(mzml_path),
            str(fasta_path),
            str(output_path),
            precursor_tol=10.0,
            fragment_tol=10.0,
        )

        assert output_path.exists()

        # Verify the file is valid idXML by loading it
        import pyopenms as oms

        prot_ids = []
        pep_ids = oms.PeptideIdentificationList()
        oms.IdXMLFile().load(str(output_path), prot_ids, pep_ids)
        assert isinstance(prot_ids, list)
