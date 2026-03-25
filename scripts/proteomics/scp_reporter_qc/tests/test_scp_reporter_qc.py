"""Tests for scp_reporter_qc."""

import csv
import math
import sys

from conftest import requires_pyopenms


@requires_pyopenms
def test_compute_sample_to_carrier_ratios():
    from scp_reporter_qc import compute_sample_to_carrier_ratios

    spectra = [
        {"spectrum_id": "s1", "126": 100.0, "127N": 120.0, "131C": 50000.0},
        {"spectrum_id": "s2", "126": 200.0, "127N": 180.0, "131C": 40000.0},
    ]
    results = compute_sample_to_carrier_ratios(spectra, "131C")
    assert len(results) == 2
    # mean of 100,120 = 110, ratio = 110/50000
    assert abs(results[0]["sample_to_carrier_ratio"] - 110.0 / 50000.0) < 1e-6
    assert results[0]["num_nonzero_samples"] == 2


@requires_pyopenms
def test_zero_carrier():
    from scp_reporter_qc import compute_sample_to_carrier_ratios

    spectra = [{"spectrum_id": "s1", "126": 100.0, "131C": 0.0}]
    results = compute_sample_to_carrier_ratios(spectra, "131C")
    assert math.isnan(results[0]["sample_to_carrier_ratio"])


@requires_pyopenms
def test_qc_summary():
    from scp_reporter_qc import qc_summary

    ratios = [
        {"sample_to_carrier_ratio": 0.002},
        {"sample_to_carrier_ratio": 0.003},
        {"sample_to_carrier_ratio": 0.005},
        {"sample_to_carrier_ratio": 0.001},
    ]
    summary = qc_summary(ratios)
    assert summary["n_spectra"] == 4
    assert abs(summary["mean_ratio"] - 0.00275) < 1e-6
    assert summary["below_0_01_count"] == 4  # all below 0.01


@requires_pyopenms
def test_qc_summary_empty():
    from scp_reporter_qc import qc_summary

    summary = qc_summary([])
    assert summary["n_spectra"] == 0
    assert math.isnan(summary["median_ratio"])


@requires_pyopenms
def test_cli_roundtrip(tmp_path):
    from scp_reporter_qc import main

    input_file = tmp_path / "input.tsv"
    output_file = tmp_path / "output.tsv"

    with open(input_file, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["spectrum_id", "126", "127N", "131C"])
        writer.writerow(["s1", "100.0", "120.0", "50000.0"])
        writer.writerow(["s2", "200.0", "180.0", "40000.0"])

    sys.argv = [
        "scp_reporter_qc.py",
        "--input", str(input_file),
        "--carrier-channel", "131C",
        "--output", str(output_file),
    ]
    main()
    assert output_file.exists()
