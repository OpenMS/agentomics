"""Tests for idxml_to_tsv_exporter."""

import os
import tempfile

from conftest import requires_pyopenms


@requires_pyopenms
def test_create_synthetic_idxml():
    import pyopenms as oms
    from idxml_to_tsv_exporter import create_synthetic_idxml

    with tempfile.TemporaryDirectory() as tmp:
        idxml_path = os.path.join(tmp, "test.idXML")
        create_synthetic_idxml(idxml_path)

        protein_ids = []
        peptide_ids = []
        oms.IdXMLFile().load(idxml_path, protein_ids, peptide_ids)
        assert len(protein_ids) == 1
        assert len(peptide_ids) == 3


@requires_pyopenms
def test_export_idxml():
    from idxml_to_tsv_exporter import create_synthetic_idxml, export_idxml

    with tempfile.TemporaryDirectory() as tmp:
        idxml_path = os.path.join(tmp, "test.idXML")
        tsv_path = os.path.join(tmp, "results.tsv")

        create_synthetic_idxml(idxml_path)
        stats = export_idxml(idxml_path, tsv_path)

        assert stats["peptide_ids"] == 3
        assert stats["total_psms"] == 3
        assert stats["protein_ids"] == 1

        with open(tsv_path) as fh:
            lines = fh.readlines()
        assert lines[0].strip().startswith("spectrum_reference")
        assert len(lines) == 4  # header + 3 PSMs


@requires_pyopenms
def test_export_content():
    from idxml_to_tsv_exporter import create_synthetic_idxml, export_idxml

    with tempfile.TemporaryDirectory() as tmp:
        idxml_path = os.path.join(tmp, "test.idXML")
        tsv_path = os.path.join(tmp, "results.tsv")

        create_synthetic_idxml(idxml_path)
        export_idxml(idxml_path, tsv_path)

        with open(tsv_path) as fh:
            lines = fh.readlines()

        # Check first data row contains expected peptide sequence
        data_row = lines[1].split("\t")
        assert len(data_row) == 8
        # sequence column is index 3
        assert data_row[3].strip() in ["ACDEFGHIK", "MNPQRSTWY"]
