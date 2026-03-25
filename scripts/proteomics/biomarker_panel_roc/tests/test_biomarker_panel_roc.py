"""Tests for biomarker_panel_roc."""

import csv
import sys

from conftest import requires_pyopenms


@requires_pyopenms
def test_compute_roc_perfect():
    from biomarker_panel_roc import compute_roc

    # Perfect separation: all cases have higher scores
    scores = [10, 9, 8, 7, 1, 2, 3, 4]
    labels = [1, 1, 1, 1, 0, 0, 0, 0]
    fpr, tpr, auc = compute_roc(scores, labels)
    assert abs(auc - 1.0) < 0.01


@requires_pyopenms
def test_compute_roc_random():
    from biomarker_panel_roc import compute_roc

    # No separation
    scores = [1, 2, 3, 4, 5, 6, 7, 8]
    labels = [1, 0, 1, 0, 1, 0, 1, 0]
    _, _, auc = compute_roc(scores, labels)
    assert 0.3 < auc < 0.8  # near 0.5


@requires_pyopenms
def test_compute_roc_all_same_class():
    from biomarker_panel_roc import compute_roc

    scores = [1, 2, 3]
    labels = [1, 1, 1]
    _, _, auc = compute_roc(scores, labels)
    assert auc == 0.5  # fallback


@requires_pyopenms
def test_analyze_biomarkers():
    from biomarker_panel_roc import analyze_biomarkers

    quant = {
        "P1": {"case_1": 100, "case_2": 110, "control_1": 50, "control_2": 45},
        "P2": {"case_1": 10, "case_2": 12, "control_1": 10, "control_2": 11},
    }
    groups = {"case_1": 1, "case_2": 1, "control_1": 0, "control_2": 0}
    results = analyze_biomarkers(quant, groups)
    assert len(results) == 2
    # P1 should have higher AUC (better separation)
    assert results[0]["protein_id"] == "P1"
    assert results[0]["auc"] > results[1]["auc"]


@requires_pyopenms
def test_cli_roundtrip(tmp_path):
    from biomarker_panel_roc import main

    input_file = tmp_path / "input.tsv"
    output_file = tmp_path / "output.tsv"

    with open(input_file, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["protein_id", "case_1", "case_2", "control_1", "control_2"])
        writer.writerow(["P1", "100", "110", "50", "45"])
        writer.writerow(["P2", "10", "12", "10", "11"])

    sys.argv = [
        "biomarker_panel_roc.py",
        "--input", str(input_file),
        "--groups", "case,control",
        "--output", str(output_file),
    ]
    main()
    assert output_file.exists()
