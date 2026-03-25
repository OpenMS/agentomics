"""Tests for immunopeptidome_qc."""

import csv

import pytest

pytest.importorskip("pyopenms")


def test_validate_sequence():
    from immunopeptidome_qc import validate_sequence

    assert validate_sequence("AAGIGILTV") is True
    assert validate_sequence("") is False


def test_length_distribution():
    from immunopeptidome_qc import length_distribution

    seqs = ["AAGIGILTV", "GILGFVFTL", "PEPTIDEK", "AAFGIILPKQR"]
    dist = length_distribution(seqs)
    assert dist[9] == 2  # two 9-mers
    assert dist[8] == 1  # one 8-mer
    assert dist[11] == 1  # one 11-mer


def test_length_qc_class_i():
    from immunopeptidome_qc import length_qc

    dist = {8: 5, 9: 50, 10: 30, 11: 10, 5: 3, 15: 2}
    qc = length_qc(dist, "I")
    assert qc["expected_min"] == 8
    assert qc["expected_max"] == 12
    assert qc["in_range_count"] == 95  # 5+50+30+10
    assert qc["out_of_range_count"] == 5  # 3+2
    assert abs(qc["in_range_fraction"] - 0.95) < 0.01


def test_length_qc_class_ii():
    from immunopeptidome_qc import length_qc

    dist = {13: 20, 15: 30, 9: 5}
    qc = length_qc(dist, "II")
    assert qc["expected_min"] == 12
    assert qc["expected_max"] == 25
    assert qc["in_range_count"] == 50
    assert qc["out_of_range_count"] == 5


def test_positional_frequencies():
    from immunopeptidome_qc import positional_frequencies

    # All 9-mers starting with A
    seqs = ["AAGIGILTV", "AAGIGILTV", "AAGIGILTV"]
    freqs = positional_frequencies(seqs, 9)
    assert len(freqs) == 9
    assert freqs[0]["A"] == 1.0  # all start with A


def test_information_content():
    from immunopeptidome_qc import information_content

    # Uniform distribution = 0 information content
    uniform = {chr(65 + i): 0.05 for i in range(20)}
    ic = information_content(uniform)
    assert abs(ic) < 0.01

    # Single residue = max information content
    single = {"A": 1.0}
    ic_max = information_content(single)
    assert ic_max > 4.0  # log2(20) ~ 4.32


def test_run_qc():
    from immunopeptidome_qc import run_qc

    seqs = ["AAGIGILTV", "GILGFVFTL", "PEPTIDEK", "AAFGIILPK"]
    dist, qc, anchors, ic_values = run_qc(seqs, "I")
    assert sum(dist.values()) == 4
    assert qc["in_range_fraction"] == 1.0  # all are 8-9aa
    assert len(ic_values) > 0


def test_cli_roundtrip(tmp_path):
    from click.testing import CliRunner
    from immunopeptidome_qc import main

    input_file = tmp_path / "input.tsv"
    output_file = tmp_path / "length_dist.tsv"
    motifs_file = tmp_path / "anchor_freq.tsv"

    with open(input_file, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["sequence"])
        for seq in ["AAGIGILTV", "GILGFVFTL", "PEPTIDEK", "AAFGIILPK"]:
            writer.writerow([seq])

    runner = CliRunner()
    result = runner.invoke(main, [
        "--input", str(input_file),
        "--hla-class", "I",
        "--output", str(output_file),
        "--motifs", str(motifs_file),
    ])
    assert result.exit_code == 0, f"CLI failed: {result.output}\n{result.exception}"

    assert output_file.exists()
    assert motifs_file.exists()
