"""Tests for metapeptide_function_aggregator."""

import csv
import sys

from conftest import requires_pyopenms


@requires_pyopenms
def test_load_peptide_protein_map(tmp_path):
    from metapeptide_function_aggregator import load_peptide_protein_map

    pep_file = tmp_path / "peptides.tsv"
    with open(pep_file, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["peptide", "protein"])
        writer.writerow(["PEPTIDEK", "P1"])
        writer.writerow(["PEPTIDEK", "P2"])
        writer.writerow(["AGIILTK", "P3"])

    pep_map = load_peptide_protein_map(str(pep_file))
    assert "PEPTIDEK" in pep_map
    assert pep_map["PEPTIDEK"] == {"P1", "P2"}
    assert pep_map["AGIILTK"] == {"P3"}


@requires_pyopenms
def test_load_peptide_protein_map_semicolon(tmp_path):
    from metapeptide_function_aggregator import load_peptide_protein_map

    pep_file = tmp_path / "peptides.tsv"
    with open(pep_file, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["peptide", "protein"])
        writer.writerow(["PEPTIDEK", "P1;P2"])

    pep_map = load_peptide_protein_map(str(pep_file))
    assert pep_map["PEPTIDEK"] == {"P1", "P2"}


@requires_pyopenms
def test_load_annotations(tmp_path):
    from metapeptide_function_aggregator import load_annotations

    ann_file = tmp_path / "annotations.tsv"
    with open(ann_file, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["protein", "term_id", "term_name"])
        writer.writerow(["P1", "GO:0006412", "translation"])
        writer.writerow(["P1", "GO:0005840", "ribosome"])

    annotations = load_annotations(str(ann_file))
    assert "P1" in annotations
    assert len(annotations["P1"]) == 2


@requires_pyopenms
def test_aggregate_functions():
    from metapeptide_function_aggregator import aggregate_functions

    pep_to_prot = {
        "PEPTIDEK": {"P1", "P2"},
        "AGIILTK": {"P3"},
    }
    prot_to_terms = {
        "P1": [("GO:0006412", "translation")],
        "P2": [("GO:0006412", "translation"), ("GO:0005840", "ribosome")],
        "P3": [("KEGG:00010", "glycolysis")],
    }

    pep_annots, term_counts = aggregate_functions(pep_to_prot, prot_to_terms)
    assert len(pep_annots) == 2

    # PEPTIDEK maps to P1 and P2, should get 2 unique terms
    peptidek_entry = [pa for pa in pep_annots if pa["peptide"] == "PEPTIDEK"][0]
    assert peptidek_entry["n_terms"] == 2

    assert term_counts["GO:0006412"] == 1  # appears once for PEPTIDEK (deduplicated)
    assert term_counts["KEGG:00010"] == 1


@requires_pyopenms
def test_summarize_terms():
    from collections import Counter

    from metapeptide_function_aggregator import summarize_terms

    term_counts = Counter({"GO:0006412": 5, "GO:0005840": 3})
    prot_to_terms = {
        "P1": [("GO:0006412", "translation"), ("GO:0005840", "ribosome")],
    }

    summary = summarize_terms(term_counts, prot_to_terms)
    assert len(summary) == 2
    assert summary[0]["term_id"] == "GO:0006412"
    assert summary[0]["peptide_count"] == 5


@requires_pyopenms
def test_cli_roundtrip(tmp_path):
    from metapeptide_function_aggregator import main

    pep_file = tmp_path / "peptides.tsv"
    ann_file = tmp_path / "annotations.tsv"
    output_file = tmp_path / "output.tsv"

    with open(pep_file, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["peptide", "protein"])
        writer.writerow(["PEPTIDEK", "P1"])
        writer.writerow(["AGIILTK", "P2"])

    with open(ann_file, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["protein", "term_id", "term_name"])
        writer.writerow(["P1", "GO:0006412", "translation"])
        writer.writerow(["P2", "KEGG:00010", "glycolysis"])

    sys.argv = [
        "metapeptide_function_aggregator.py",
        "--peptides", str(pep_file),
        "--annotations", str(ann_file),
        "--output", str(output_file),
    ]
    main()
    assert output_file.exists()
