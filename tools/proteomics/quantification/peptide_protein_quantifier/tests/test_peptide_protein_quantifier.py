"""Tests for peptide_protein_quantifier."""

import csv
import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestCreateSyntheticAnnotatedConsensus:
    def test_creates_valid_file(self):
        import pyopenms as oms
        from peptide_protein_quantifier import create_synthetic_annotated_consensus

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "consensus.consensusXML")
            create_synthetic_annotated_consensus(path)

            cm = oms.ConsensusMap()
            oms.ConsensusXMLFile().load(path, cm)
            assert cm.size() == 4  # 2 proteins x 2 peptides

    def test_has_peptide_annotations(self):
        import pyopenms as oms
        from peptide_protein_quantifier import create_synthetic_annotated_consensus

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "consensus.consensusXML")
            create_synthetic_annotated_consensus(path)

            cm = oms.ConsensusMap()
            oms.ConsensusXMLFile().load(path, cm)

            annotated = 0
            for cf in cm:
                pids = cf.getPeptideIdentifications()
                if pids.size() > 0:
                    annotated += 1
            assert annotated == 4

    def test_has_protein_identifications(self):
        import pyopenms as oms
        from peptide_protein_quantifier import create_synthetic_annotated_consensus

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "consensus.consensusXML")
            create_synthetic_annotated_consensus(path)

            cm = oms.ConsensusMap()
            oms.ConsensusXMLFile().load(path, cm)

            prot_ids = cm.getProteinIdentifications()
            assert len(prot_ids) == 1
            assert len(prot_ids[0].getHits()) == 2

    def test_custom_proteins(self):
        import pyopenms as oms
        from peptide_protein_quantifier import create_synthetic_annotated_consensus

        custom = {
            "MyProtein": [
                ("TESTPEPTIDEK", 450.0, 50.0, 5000.0),
            ],
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "consensus.consensusXML")
            create_synthetic_annotated_consensus(path, proteins=custom)

            cm = oms.ConsensusMap()
            oms.ConsensusXMLFile().load(path, cm)
            assert cm.size() == 1


class TestQuantifyProteins:
    def test_two_protein_rollup(self):
        from peptide_protein_quantifier import (
            create_synthetic_annotated_consensus,
            quantify_proteins,
        )

        with tempfile.TemporaryDirectory() as tmp:
            in_path = os.path.join(tmp, "input.consensusXML")
            out_path = os.path.join(tmp, "protein_quant.csv")

            create_synthetic_annotated_consensus(in_path)

            n_proteins = quantify_proteins(in_path, out_path, top_n=3, include_all=True)
            assert n_proteins == 2
            assert os.path.exists(out_path)

    def test_output_csv_content(self):
        from peptide_protein_quantifier import (
            create_synthetic_annotated_consensus,
            quantify_proteins,
        )

        with tempfile.TemporaryDirectory() as tmp:
            in_path = os.path.join(tmp, "input.consensusXML")
            out_path = os.path.join(tmp, "protein_quant.csv")

            create_synthetic_annotated_consensus(in_path)
            quantify_proteins(in_path, out_path)

            with open(out_path) as f:
                reader = csv.reader(f)
                rows = list(reader)

            # Header row
            assert rows[0] == ["protein_accession", "score"]
            # Should have protein rows
            accessions = {row[0] for row in rows[1:] if row and not row[0].startswith("#")}
            assert "ProteinA" in accessions
            assert "ProteinB" in accessions

    def test_single_protein(self):
        from peptide_protein_quantifier import (
            create_synthetic_annotated_consensus,
            quantify_proteins,
        )

        custom = {
            "SingleProt": [
                ("PEPK", 500.0, 100.0, 10000.0),
                ("SEQR", 600.0, 200.0, 20000.0),
                ("ANOTHERK", 700.0, 300.0, 30000.0),
            ],
        }
        with tempfile.TemporaryDirectory() as tmp:
            in_path = os.path.join(tmp, "input.consensusXML")
            out_path = os.path.join(tmp, "protein_quant.csv")

            create_synthetic_annotated_consensus(in_path, proteins=custom)
            n_proteins = quantify_proteins(in_path, out_path, top_n=3, include_all=True)
            assert n_proteins == 1
