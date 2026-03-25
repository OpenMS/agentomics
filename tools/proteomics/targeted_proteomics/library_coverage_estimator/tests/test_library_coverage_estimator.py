"""Tests for library_coverage_estimator."""

import csv

from conftest import requires_pyopenms


@requires_pyopenms
def test_read_library_peptides(tmp_path):
    from library_coverage_estimator import read_library_peptides

    lib_file = tmp_path / "lib.tsv"
    with open(lib_file, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["PeptideSequence", "charge"])
        writer.writerow(["PEPTIDEK", "2"])
        writer.writerow(["AGIILTK", "2"])
        writer.writerow(["PEPTIDEK", "3"])  # duplicate sequence

    peps = read_library_peptides(str(lib_file))
    assert "PEPTIDEK" in peps
    assert "AGIILTK" in peps
    assert len(peps) == 2  # deduplicated


@requires_pyopenms
def test_digest_fasta(tmp_path):
    import pyopenms as oms
    from library_coverage_estimator import digest_fasta

    fasta_file = tmp_path / "test.fasta"
    entries = [oms.FASTAEntry()]
    entries[0].identifier = "P12345"
    entries[0].sequence = "MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVK"
    oms.FASTAFile().store(str(fasta_file), entries)

    prot_peps, all_peps = digest_fasta(str(fasta_file), "Trypsin", 1)
    assert "P12345" in prot_peps
    assert len(all_peps) > 0


@requires_pyopenms
def test_compute_coverage():
    from library_coverage_estimator import compute_coverage

    lib_peps = {"PEPTIDEK", "AGIILTK", "UNKNOWN"}
    prot_peps = {
        "P1": ["PEPTIDEK", "AGIILTK", "FOOBAR"],
        "P2": ["XYZTHING"],
    }
    all_digest = {"PEPTIDEK", "AGIILTK", "FOOBAR", "XYZTHING"}

    result = compute_coverage(lib_peps, prot_peps, all_digest)
    assert result["total_digestible_peptides"] == 4
    assert result["library_peptides_in_proteome"] == 2  # PEPTIDEK, AGIILTK
    assert result["proteins_with_library_peptide"] == 1  # only P1
    assert result["total_proteins"] == 2


@requires_pyopenms
def test_cli_roundtrip(tmp_path):
    import pyopenms as oms
    from click.testing import CliRunner
    from library_coverage_estimator import main

    # Create FASTA
    fasta_file = tmp_path / "proteome.fasta"
    entries = [oms.FASTAEntry()]
    entries[0].identifier = "P12345"
    entries[0].sequence = "MKWVTFISLLFLFSSAYSRGVFRRDAHK"
    oms.FASTAFile().store(str(fasta_file), entries)

    # Digest to find real peptides
    from library_coverage_estimator import digest_fasta
    prot_peps, all_peps = digest_fasta(str(fasta_file), "Trypsin", 1)
    some_peps = list(all_peps)[:2] if all_peps else ["PEPTIDEK"]

    # Create library
    lib_file = tmp_path / "lib.tsv"
    with open(lib_file, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["PeptideSequence"])
        for p in some_peps:
            writer.writerow([p])

    output_file = tmp_path / "coverage.tsv"
    runner = CliRunner()
    result = runner.invoke(main, [
        "--library", str(lib_file),
        "--fasta", str(fasta_file),
        "--enzyme", "Trypsin",
        "--output", str(output_file),
    ])
    assert result.exit_code == 0, f"CLI failed: {result.output}\n{result.exception}"
    assert output_file.exists()
