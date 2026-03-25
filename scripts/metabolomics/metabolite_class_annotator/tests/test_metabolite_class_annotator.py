"""Tests for metabolite_class_annotator."""

from conftest import requires_pyopenms


@requires_pyopenms
class TestMetaboliteClassAnnotator:
    def test_compute_mass_defect(self):
        from metabolite_class_annotator import compute_mass_defect

        md = compute_mass_defect(180.0634)
        assert abs(md - 0.0634) < 0.001

    def test_compute_kendrick_md(self):
        from metabolite_class_annotator import compute_kendrick_mass_defect

        kmd = compute_kendrick_mass_defect(180.0634)
        assert isinstance(kmd, float)

    def test_annotate_class_small_molecule(self):
        from metabolite_class_annotator import annotate_class

        classes = annotate_class(180.0634)  # glucose
        assert len(classes) > 0
        assert isinstance(classes[0], str)

    def test_annotate_class_lipid_range(self):
        from metabolite_class_annotator import annotate_class

        classes = annotate_class(700.2)
        assert "Lipid" in classes

    def test_annotate_features(self):
        from metabolite_class_annotator import annotate_features

        features = [
            {"mz": "181.0707", "rt": "60.0", "intensity": "1000"},
            {"mz": "701.2", "rt": "300.0", "intensity": "5000"},
        ]
        results = annotate_features(features)
        assert len(results) == 2
        assert "compound_class" in results[0]
        assert "mass_defect" in results[0]
        assert "kendrick_md" in results[0]

    def test_unknown_class(self):
        from metabolite_class_annotator import annotate_class

        # Very large mass outside normal ranges
        classes = annotate_class(50000.0)
        assert "Unknown" in classes
