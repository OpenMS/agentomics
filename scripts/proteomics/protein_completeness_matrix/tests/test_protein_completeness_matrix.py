"""Tests for protein_completeness_matrix."""

import csv
import sys

import numpy as np
from conftest import requires_pyopenms


@requires_pyopenms
def test_load_quant_matrix(tmp_path):
    from protein_completeness_matrix import load_quant_matrix

    input_file = tmp_path / "quant.tsv"
    with open(input_file, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["protein_id", "s1", "s2", "s3"])
        writer.writerow(["P1", "100.0", "NA", "95.0"])
        writer.writerow(["P2", "200.0", "180.0", "190.0"])

    pids, sids, data = load_quant_matrix(str(input_file))
    assert pids == ["P1", "P2"]
    assert sids == ["s1", "s2", "s3"]
    assert data.shape == (2, 3)
    assert np.isnan(data[0, 1])  # P1, s2 is NA
    assert data[1, 0] == 200.0


@requires_pyopenms
def test_compute_protein_completeness():
    from protein_completeness_matrix import compute_protein_completeness

    data = np.array([
        [1.0, np.nan, 1.0],  # 2/3 complete
        [1.0, 1.0, 1.0],    # 3/3 complete
    ])
    comp = compute_protein_completeness(data)
    assert abs(comp[0] - 2.0 / 3.0) < 0.01
    assert abs(comp[1] - 1.0) < 0.01


@requires_pyopenms
def test_compute_sample_completeness():
    from protein_completeness_matrix import compute_sample_completeness

    data = np.array([
        [1.0, np.nan, 1.0],
        [1.0, 1.0, 1.0],
    ])
    comp = compute_sample_completeness(data)
    assert abs(comp[0] - 1.0) < 0.01     # s1: both proteins present
    assert abs(comp[1] - 0.5) < 0.01     # s2: one of two
    assert abs(comp[2] - 1.0) < 0.01     # s3: both present


@requires_pyopenms
def test_filter_by_completeness():
    from protein_completeness_matrix import filter_by_completeness

    data = np.array([
        [1.0, np.nan, np.nan],  # 1/3 = 0.33
        [1.0, 1.0, 1.0],       # 3/3 = 1.0
        [1.0, 1.0, np.nan],    # 2/3 = 0.67
    ])
    comp = np.array([1.0 / 3, 1.0, 2.0 / 3])
    filtered_ids, filtered_data = filter_by_completeness(
        ["P1", "P2", "P3"], data, comp, 0.5
    )
    assert filtered_ids == ["P2", "P3"]
    assert filtered_data.shape == (2, 3)


@requires_pyopenms
def test_completeness_summary():
    from protein_completeness_matrix import completeness_summary

    data = np.array([
        [1.0, np.nan],
        [1.0, 1.0],
    ])
    prot_comp = np.array([0.5, 1.0])
    samp_comp = np.array([1.0, 0.5])
    summary = completeness_summary(prot_comp, samp_comp, data)
    assert summary["total_proteins"] == 2
    assert summary["total_samples"] == 2
    assert summary["non_missing_cells"] == 3
    assert abs(summary["overall_completeness"] - 0.75) < 0.01


@requires_pyopenms
def test_cli_roundtrip(tmp_path):
    from protein_completeness_matrix import main

    input_file = tmp_path / "quant.tsv"
    output_file = tmp_path / "completeness.tsv"

    with open(input_file, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["protein_id", "s1", "s2", "s3"])
        writer.writerow(["P1", "100.0", "NA", "95.0"])
        writer.writerow(["P2", "200.0", "180.0", "190.0"])
        writer.writerow(["P3", "NA", "NA", "50.0"])

    sys.argv = [
        "protein_completeness_matrix.py",
        "--input", str(input_file),
        "--min-completeness", "0.5",
        "--output", str(output_file),
    ]
    main()
    assert output_file.exists()
