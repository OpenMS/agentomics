"""Tests for rt_alignment_identification."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestAlignIdentification:
    def test_align_removes_rt_offset(self):
        """Two idXML files with +10.0 RT offset should align to < 1.0 difference."""
        import pyopenms as oms
        from rt_alignment_identification import (
            align_identification,
            create_synthetic_idxml,
        )

        with tempfile.TemporaryDirectory() as tmp:
            ref_path = os.path.join(tmp, "reference.idXML")
            inp_path = os.path.join(tmp, "input.idXML")
            out_path = os.path.join(tmp, "aligned.idXML")

            create_synthetic_idxml(ref_path, n_ids=5, rt_offset=0.0)
            create_synthetic_idxml(inp_path, n_ids=5, rt_offset=10.0)

            stats = align_identification(ref_path, inp_path, out_path, model="linear")

            assert os.path.exists(out_path)
            assert stats["n_peptide_ids"] == 5

            # Load reference and aligned, compare RTs
            ref_prot = []
            ref_pep = oms.PeptideIdentificationList()
            oms.IdXMLFile().load(ref_path, ref_prot, ref_pep)

            aligned_prot = []
            aligned_pep = oms.PeptideIdentificationList()
            oms.IdXMLFile().load(out_path, aligned_prot, aligned_pep)

            for i in range(ref_pep.size()):
                rt_diff = abs(ref_pep.at(i).getRT() - aligned_pep.at(i).getRT())
                assert rt_diff < 1.0, (
                    f"Peptide ID {i}: RT diff {rt_diff:.2f} >= 1.0"
                )

    def test_align_returns_stats(self):
        """Stats dict should contain expected keys."""
        from rt_alignment_identification import (
            align_identification,
            create_synthetic_idxml,
        )

        with tempfile.TemporaryDirectory() as tmp:
            ref_path = os.path.join(tmp, "reference.idXML")
            inp_path = os.path.join(tmp, "input.idXML")
            out_path = os.path.join(tmp, "aligned.idXML")

            create_synthetic_idxml(ref_path, n_ids=5, rt_offset=0.0)
            create_synthetic_idxml(inp_path, n_ids=5, rt_offset=5.0)

            stats = align_identification(ref_path, inp_path, out_path)

            assert "n_peptide_ids" in stats
            assert "n_trafo_points" in stats
            assert "model" in stats
            assert isinstance(stats["n_peptide_ids"], int)

    def test_align_preserves_peptide_count(self):
        """Alignment should not change the number of peptide identifications."""
        import pyopenms as oms
        from rt_alignment_identification import (
            align_identification,
            create_synthetic_idxml,
        )

        with tempfile.TemporaryDirectory() as tmp:
            ref_path = os.path.join(tmp, "reference.idXML")
            inp_path = os.path.join(tmp, "input.idXML")
            out_path = os.path.join(tmp, "aligned.idXML")

            create_synthetic_idxml(ref_path, n_ids=5, rt_offset=0.0)
            create_synthetic_idxml(inp_path, n_ids=5, rt_offset=15.0)

            align_identification(ref_path, inp_path, out_path)

            aligned_prot = []
            aligned_pep = oms.PeptideIdentificationList()
            oms.IdXMLFile().load(out_path, aligned_prot, aligned_pep)

            assert aligned_pep.size() == 5

    def test_create_synthetic_idxml(self):
        """Synthetic idXML should have correct number of peptide IDs."""
        import pyopenms as oms
        from rt_alignment_identification import create_synthetic_idxml

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "test.idXML")
            create_synthetic_idxml(path, n_ids=4, rt_offset=20.0)

            prot_ids = []
            pep_ids = oms.PeptideIdentificationList()
            oms.IdXMLFile().load(path, prot_ids, pep_ids)

            assert pep_ids.size() == 4
            assert pep_ids.at(0).getRT() == pytest.approx(120.0, abs=0.01)
