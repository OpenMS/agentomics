# AI Contributor Skills & Validation Pipeline — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create skills, contributor docs, CI, and per-script directory structure so AI agents can contribute validated pyopenms scripts to agentomics.

**Architecture:** Restructure existing flat scripts into self-contained per-script directories with isolated dependency resolution. Add Claude Code skills for guided contribution and validation, a platform-agnostic AGENTS.md, and a GitHub Actions CI pipeline that validates each script in its own venv.

**Tech Stack:** Python, pyopenms, pytest, ruff, GitHub Actions

**Spec:** `docs/superpowers/specs/2026-03-24-ai-contributor-skills-design.md`

---

### Task 1: Create ruff.toml

**Files:**
- Create: `ruff.toml`

- [ ] **Step 1: Create ruff.toml**

```toml
line-length = 120

[lint]
select = ["E", "F", "W", "I"]
```

- [ ] **Step 2: Verify ruff.toml is valid**

Run: `pip install ruff && ruff check --config ruff.toml . 2>&1 | head -5`
Expected: No config errors (may show "No files found" which is fine)

- [ ] **Step 3: Commit**

```bash
git add ruff.toml
git commit -m "Add ruff.toml with E/F/W/I rule set, line-length 120"
```

---

### Task 2: Migrate peptide_mass_calculator to per-script directory

**Files:**
- Create: `tools/proteomics/peptide_mass_calculator/peptide_mass_calculator.py`
- Create: `tools/proteomics/peptide_mass_calculator/requirements.txt`
- Create: `tools/proteomics/peptide_mass_calculator/README.md`
- Create: `tools/proteomics/peptide_mass_calculator/tests/conftest.py`
- Create: `tools/proteomics/peptide_mass_calculator/tests/test_peptide_mass_calculator.py`

- [ ] **Step 1: Create directory structure**

```bash
mkdir -p tools/proteomics/peptide_mass_calculator/tests
```

- [ ] **Step 2: Copy script from feature branch**

```bash
git show origin/copilot/add-agentic-scripts-for-proteomics:tools/proteomics/peptide_mass_calculator.py > tools/proteomics/peptide_mass_calculator/peptide_mass_calculator.py
```

- [ ] **Step 3: Create requirements.txt**

```
pyopenms
```

- [ ] **Step 4: Create conftest.py**

```python
import sys
import os

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

try:
    import pyopenms  # noqa: F401

    HAS_PYOPENMS = True
except ImportError:
    HAS_PYOPENMS = False

requires_pyopenms = pytest.mark.skipif(not HAS_PYOPENMS, reason="pyopenms not installed")
```

- [ ] **Step 5: Create test file**

Extract the `TestPeptideMassCalculator` class from the feature branch's `tests/test_proteomics.py`. Adapt imports to use the conftest marker and direct module imports (no sys.path manipulation in the test file itself):

```python
"""Tests for peptide_mass_calculator."""

import pytest

from conftest import requires_pyopenms


@requires_pyopenms
class TestPeptideMassCalculator:
    def test_basic_mass(self):
        from peptide_mass_calculator import peptide_masses

        result = peptide_masses("PEPTIDEK")
        assert result["sequence"] == "PEPTIDEK"
        assert result["charge"] == 1
        assert 927.0 < result["monoisotopic_mass"] < 928.0
        assert result["mz_monoisotopic"] > result["monoisotopic_mass"]

    def test_charge_state(self):
        from peptide_mass_calculator import peptide_masses

        r1 = peptide_masses("PEPTIDEK", charge=1)
        r2 = peptide_masses("PEPTIDEK", charge=2)
        assert r2["mz_monoisotopic"] < r1["mz_monoisotopic"]

    def test_fragment_ions(self):
        from peptide_mass_calculator import fragment_ions

        ions = fragment_ions("PEPTIDEK")
        seq_len = len("PEPTIDEK")
        assert len(ions["b_ions"]) == seq_len - 1
        assert len(ions["y_ions"]) == seq_len - 1

    def test_modified_sequence(self):
        from peptide_mass_calculator import peptide_masses

        result = peptide_masses("PEPTM[147]IDEK")
        assert result["monoisotopic_mass"] > 0

    def test_mz_formula(self):
        from peptide_mass_calculator import peptide_masses, PROTON

        r = peptide_masses("PEPTIDEK", charge=2)
        expected = (r["monoisotopic_mass"] + 2 * PROTON) / 2
        assert abs(r["mz_monoisotopic"] - expected) < 1e-6
```

- [ ] **Step 6: Create README.md**

```markdown
# Peptide Mass Calculator

Calculate monoisotopic and average masses for peptide sequences, and compute
b-ion / y-ion fragment series.

## Usage

```bash
python peptide_mass_calculator.py --sequence PEPTIDEK
python peptide_mass_calculator.py --sequence PEPTM[147]IDEK --charge 2
python peptide_mass_calculator.py --sequence ACDEFGHIK --fragments
```
```

- [ ] **Step 7: Run ruff**

Run: `ruff check tools/proteomics/peptide_mass_calculator/`
Expected: No errors

- [ ] **Step 8: Run tests**

Run: `PYTHONPATH=tools/proteomics/peptide_mass_calculator python -m pytest tools/proteomics/peptide_mass_calculator/tests/ -v`
Expected: 5 tests pass (or skip if pyopenms not installed)

- [ ] **Step 9: Commit**

```bash
git add tools/proteomics/peptide_mass_calculator/
git commit -m "Migrate peptide_mass_calculator to per-script directory structure"
```

---

### Task 3: Migrate protein_digest to per-script directory

**Files:**
- Create: `tools/proteomics/protein_digest/protein_digest.py`
- Create: `tools/proteomics/protein_digest/requirements.txt`
- Create: `tools/proteomics/protein_digest/README.md`
- Create: `tools/proteomics/protein_digest/tests/conftest.py`
- Create: `tools/proteomics/protein_digest/tests/test_protein_digest.py`

- [ ] **Step 1: Create directory and copy script**

```bash
mkdir -p tools/proteomics/protein_digest/tests
git show origin/copilot/add-agentic-scripts-for-proteomics:tools/proteomics/protein_digest.py > tools/proteomics/protein_digest/protein_digest.py
```

- [ ] **Step 2: Create requirements.txt**

```
pyopenms
```

- [ ] **Step 3: Create conftest.py** (identical to Task 2 Step 4)

- [ ] **Step 4: Create test file**

Extract the `TestProteinDigest` class from the feature branch's `tests/test_proteomics.py`:

```python
"""Tests for protein_digest."""

import pytest

from conftest import requires_pyopenms


@requires_pyopenms
class TestProteinDigest:
    PROTEIN = "MKVLWAALLVTFLAGCQAKVEQAVETEPEPELRQQTEWQSGQRWELAL"

    def test_tryptic_digest_returns_peptides(self):
        from protein_digest import digest_protein

        peptides = digest_protein(self.PROTEIN, enzyme="Trypsin", min_length=1)
        assert len(peptides) > 0

    def test_peptide_structure(self):
        from protein_digest import digest_protein

        peptides = digest_protein(self.PROTEIN, enzyme="Trypsin", min_length=1)
        for pep in peptides:
            assert "sequence" in pep
            assert "monoisotopic_mass" in pep
            assert pep["monoisotopic_mass"] > 0

    def test_missed_cleavages(self):
        from protein_digest import digest_protein

        peps_0 = digest_protein(self.PROTEIN, enzyme="Trypsin", missed_cleavages=0, min_length=1)
        peps_2 = digest_protein(self.PROTEIN, enzyme="Trypsin", missed_cleavages=2, min_length=1)
        assert len(peps_2) >= len(peps_0)

    def test_length_filter(self):
        from protein_digest import digest_protein

        peptides = digest_protein(
            self.PROTEIN, enzyme="Trypsin", min_length=5, max_length=20, missed_cleavages=2
        )
        for pep in peptides:
            assert 5 <= pep["length"] <= 20

    def test_list_enzymes(self):
        from protein_digest import list_enzymes

        enzymes = list_enzymes()
        assert "Trypsin" in enzymes
        assert len(enzymes) > 5
```

- [ ] **Step 5: Create README.md**

```markdown
# Protein In-Silico Digest

Perform in-silico enzymatic digestion of a protein sequence and report
the resulting peptides with their masses.

## Usage

```bash
python protein_digest.py --sequence MKVLWAALLVTFLAGCQAK... --enzyme Trypsin
python protein_digest.py --sequence MKVLWAALLVTFLAGCQAK... --enzyme Lys-C --missed-cleavages 2
python protein_digest.py --list-enzymes
```
```

- [ ] **Step 6: Run ruff and tests**

Run: `ruff check tools/proteomics/protein_digest/ && PYTHONPATH=tools/proteomics/protein_digest python -m pytest tools/proteomics/protein_digest/tests/ -v`
Expected: Lint clean, 5 tests pass

- [ ] **Step 7: Commit**

```bash
git add tools/proteomics/protein_digest/
git commit -m "Migrate protein_digest to per-script directory structure"
```

---

### Task 4: Migrate spectrum_file_info to per-script directory

**Files:**
- Create: `tools/proteomics/spectrum_file_info/spectrum_file_info.py`
- Create: `tools/proteomics/spectrum_file_info/requirements.txt`
- Create: `tools/proteomics/spectrum_file_info/README.md`
- Create: `tools/proteomics/spectrum_file_info/tests/conftest.py`
- Create: `tools/proteomics/spectrum_file_info/tests/test_spectrum_file_info.py`

- [ ] **Step 1: Create directory and copy script**

```bash
mkdir -p tools/proteomics/spectrum_file_info/tests
git show origin/copilot/add-agentic-scripts-for-proteomics:tools/proteomics/spectrum_file_info.py > tools/proteomics/spectrum_file_info/spectrum_file_info.py
```

- [ ] **Step 2: Create requirements.txt**

```
pyopenms
```

- [ ] **Step 3: Create conftest.py** (identical to Task 2 Step 4)

- [ ] **Step 4: Create test file with synthetic data**

This script processes mzML files. Tests generate synthetic MSExperiment data using pyopenms:

```python
"""Tests for spectrum_file_info."""

import pytest

from conftest import requires_pyopenms


@requires_pyopenms
class TestSpectrumFileInfo:
    def _make_experiment(self, n_spectra=5, ms_level=1):
        """Create a synthetic MSExperiment for testing."""
        import pyopenms as oms
        import numpy as np

        exp = oms.MSExperiment()
        for i in range(n_spectra):
            spec = oms.MSSpectrum()
            spec.setMSLevel(ms_level)
            spec.setRT(60.0 * i)
            mzs = np.array([100.0 + j for j in range(10)], dtype=np.float64)
            intensities = np.array([1000.0 * (j + 1) for j in range(10)], dtype=np.float64)
            spec.set_peaks([mzs, intensities])
            exp.addSpectrum(spec)
        return exp

    def test_summarise_nonempty(self):
        from spectrum_file_info import summarise_experiment

        exp = self._make_experiment(n_spectra=3)
        summary = summarise_experiment(exp)
        assert summary["n_spectra"] == 3
        assert 1 in summary["ms_levels"]

    def test_summarise_empty(self):
        from spectrum_file_info import summarise_experiment
        import pyopenms as oms

        exp = oms.MSExperiment()
        summary = summarise_experiment(exp)
        assert summary["n_spectra"] == 0

    def test_rt_range(self):
        from spectrum_file_info import summarise_experiment

        exp = self._make_experiment(n_spectra=5)
        summary = summarise_experiment(exp)
        rt_min, rt_max = summary["rt_range_sec"]
        assert rt_min == 0.0
        assert rt_max == 240.0

    def test_mz_range(self):
        from spectrum_file_info import summarise_experiment

        exp = self._make_experiment(n_spectra=2)
        summary = summarise_experiment(exp)
        mz_min, mz_max = summary["mz_range"]
        assert mz_min == pytest.approx(100.0)
        assert mz_max == pytest.approx(109.0)
```

- [ ] **Step 5: Create README.md**

```markdown
# Mass Spectrum File Info

Summarise the contents of an mzML file: spectra counts by MS level,
retention time range, m/z range, and TIC statistics.

## Usage

```bash
python spectrum_file_info.py --input sample.mzML
python spectrum_file_info.py --input sample.mzML --tic
```
```

- [ ] **Step 6: Run ruff and tests**

Run: `ruff check tools/proteomics/spectrum_file_info/ && PYTHONPATH=tools/proteomics/spectrum_file_info python -m pytest tools/proteomics/spectrum_file_info/tests/ -v`
Expected: Lint clean, 4 tests pass

- [ ] **Step 7: Commit**

```bash
git add tools/proteomics/spectrum_file_info/
git commit -m "Migrate spectrum_file_info to per-script directory with synthetic test data"
```

---

### Task 5: Migrate feature_detection_proteomics to per-script directory

**Files:**
- Create: `tools/proteomics/feature_detection_proteomics/feature_detection_proteomics.py`
- Create: `tools/proteomics/feature_detection_proteomics/requirements.txt`
- Create: `tools/proteomics/feature_detection_proteomics/README.md`
- Create: `tools/proteomics/feature_detection_proteomics/tests/conftest.py`
- Create: `tools/proteomics/feature_detection_proteomics/tests/test_feature_detection_proteomics.py`

- [ ] **Step 1: Create directory and copy script**

```bash
mkdir -p tools/proteomics/feature_detection_proteomics/tests
git show origin/copilot/add-agentic-scripts-for-proteomics:tools/proteomics/feature_detection_proteomics.py > tools/proteomics/feature_detection_proteomics/feature_detection_proteomics.py
```

- [ ] **Step 2: Create requirements.txt**

```
pyopenms
```

- [ ] **Step 3: Create conftest.py** (identical to Task 2 Step 4)

- [ ] **Step 4: Create test file with synthetic data**

This script requires mzML file input. Tests generate a synthetic MSExperiment, write it to a temp mzML file, and run feature detection on it:

```python
"""Tests for feature_detection_proteomics."""

import os
import tempfile

import pytest

from conftest import requires_pyopenms


@requires_pyopenms
class TestFeatureDetectionProteomics:
    def test_detect_features_returns_feature_map(self):
        import pyopenms as oms
        import numpy as np
        from feature_detection_proteomics import detect_features

        # Create a minimal synthetic experiment with a few peaks
        exp = oms.MSExperiment()
        for i in range(10):
            spec = oms.MSSpectrum()
            spec.setMSLevel(1)
            spec.setRT(60.0 + i * 2.0)
            mzs = np.array([500.0, 500.5, 501.0], dtype=np.float64)
            intensities = np.array([1e4, 5e3, 1e3], dtype=np.float64)
            spec.set_peaks([mzs, intensities])
            exp.addSpectrum(spec)

        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "test.mzML")
            output_path = os.path.join(tmpdir, "test.featureXML")
            oms.MzMLFile().store(input_path, exp)

            fm = detect_features(input_path, output_path)
            assert isinstance(fm, oms.FeatureMap)
            assert os.path.exists(output_path)
```

- [ ] **Step 5: Create README.md**

```markdown
# Feature Detection for Proteomics

Detect peptide isotope features in centroided LC-MS/MS data using the
`FeatureFinderCentroided` algorithm. Output is written as a featureXML file.

## Usage

```bash
python feature_detection_proteomics.py --input sample.mzML
python feature_detection_proteomics.py --input sample.mzML --output features.featureXML
```
```

- [ ] **Step 6: Run ruff and tests**

Run: `ruff check tools/proteomics/feature_detection_proteomics/ && PYTHONPATH=tools/proteomics/feature_detection_proteomics python -m pytest tools/proteomics/feature_detection_proteomics/tests/ -v`
Expected: Lint clean, 1 test passes

- [ ] **Step 7: Commit**

```bash
git add tools/proteomics/feature_detection_proteomics/
git commit -m "Migrate feature_detection_proteomics to per-script directory with synthetic test data"
```

---

### Task 6: Migrate mass_accuracy_calculator to per-script directory

**Files:**
- Create: `tools/metabolomics/mass_accuracy_calculator/mass_accuracy_calculator.py`
- Create: `tools/metabolomics/mass_accuracy_calculator/requirements.txt`
- Create: `tools/metabolomics/mass_accuracy_calculator/README.md`
- Create: `tools/metabolomics/mass_accuracy_calculator/tests/conftest.py`
- Create: `tools/metabolomics/mass_accuracy_calculator/tests/test_mass_accuracy_calculator.py`

- [ ] **Step 1: Create directory and copy script**

```bash
mkdir -p tools/metabolomics/mass_accuracy_calculator/tests
git show origin/copilot/add-agentic-scripts-for-proteomics:tools/metabolomics/mass_accuracy_calculator.py > tools/metabolomics/mass_accuracy_calculator/mass_accuracy_calculator.py
```

- [ ] **Step 2: Create requirements.txt**

```
pyopenms
```

- [ ] **Step 3: Create conftest.py** (identical to Task 2 Step 4)

- [ ] **Step 4: Create test file**

Extract `TestMassAccuracyCalculator` from the feature branch's `tests/test_metabolomics.py`:

```python
"""Tests for mass_accuracy_calculator."""

import pytest

from conftest import requires_pyopenms


@requires_pyopenms
class TestMassAccuracyCalculator:
    def test_sequence_theoretical(self):
        from mass_accuracy_calculator import theoretical_mz_from_sequence

        mz = theoretical_mz_from_sequence("PEPTIDEK", 1)
        assert 928.0 < mz < 929.0

    def test_formula_theoretical(self):
        from mass_accuracy_calculator import theoretical_mz_from_formula

        mz = theoretical_mz_from_formula("C6H12O6", 1)
        assert 181.0 < mz < 182.0

    def test_ppm_zero_error(self):
        from mass_accuracy_calculator import ppm_error

        assert ppm_error(500.0, 500.0) == 0.0

    def test_ppm_positive_error(self):
        from mass_accuracy_calculator import ppm_error

        assert ppm_error(500.0, 500.001) > 0

    def test_ppm_negative_error(self):
        from mass_accuracy_calculator import ppm_error

        assert ppm_error(500.0, 499.999) < 0

    def test_ppm_known_value(self):
        from mass_accuracy_calculator import ppm_error

        ppm = ppm_error(1000.0, 1000.001)
        assert abs(ppm - 1.0) < 0.001
```

- [ ] **Step 5: Create README.md**

```markdown
# Mass Accuracy Calculator

Calculate mass accuracy (ppm error) between a theoretical value derived
from a peptide sequence or molecular formula and observed m/z values.

## Usage

```bash
python mass_accuracy_calculator.py --sequence PEPTIDEK --observed 803.4560
python mass_accuracy_calculator.py --formula C6H12O6 --observed 181.0709
python mass_accuracy_calculator.py --sequence ACDEFGHIK --charge 2 --observed 554.2478 554.2480
```
```

- [ ] **Step 6: Run ruff and tests**

Run: `ruff check tools/metabolomics/mass_accuracy_calculator/ && PYTHONPATH=tools/metabolomics/mass_accuracy_calculator python -m pytest tools/metabolomics/mass_accuracy_calculator/tests/ -v`
Expected: Lint clean, 6 tests pass

- [ ] **Step 7: Commit**

```bash
git add tools/metabolomics/mass_accuracy_calculator/
git commit -m "Migrate mass_accuracy_calculator to per-script directory structure"
```

---

### Task 7: Migrate isotope_pattern_matcher to per-script directory

**Files:**
- Create: `tools/metabolomics/isotope_pattern_matcher/isotope_pattern_matcher.py`
- Create: `tools/metabolomics/isotope_pattern_matcher/requirements.txt`
- Create: `tools/metabolomics/isotope_pattern_matcher/README.md`
- Create: `tools/metabolomics/isotope_pattern_matcher/tests/conftest.py`
- Create: `tools/metabolomics/isotope_pattern_matcher/tests/test_isotope_pattern_matcher.py`

- [ ] **Step 1: Create directory and copy script**

```bash
mkdir -p tools/metabolomics/isotope_pattern_matcher/tests
git show origin/copilot/add-agentic-scripts-for-proteomics:tools/metabolomics/isotope_pattern_matcher.py > tools/metabolomics/isotope_pattern_matcher/isotope_pattern_matcher.py
```

- [ ] **Step 2: Create requirements.txt**

```
pyopenms
```

- [ ] **Step 3: Create conftest.py** (identical to Task 2 Step 4)

- [ ] **Step 4: Create test file**

Extract `TestIsotopePatternMatcher` from the feature branch's `tests/test_metabolomics.py`:

```python
"""Tests for isotope_pattern_matcher."""

import pytest

from conftest import requires_pyopenms


@requires_pyopenms
class TestIsotopePatternMatcher:
    def test_glucose_pattern(self):
        from isotope_pattern_matcher import get_isotope_distribution

        dist = get_isotope_distribution("C6H12O6", max_isotopes=3)
        assert len(dist) == 3
        assert dist[0][1] == pytest.approx(100.0)
        assert dist[1][1] < dist[0][1]

    def test_pattern_mz_ordering(self):
        from isotope_pattern_matcher import get_isotope_distribution

        dist = get_isotope_distribution("C12H22O11", max_isotopes=4)
        mzs = [mz for mz, _ in dist]
        assert mzs == sorted(mzs)

    def test_cosine_similarity_perfect(self):
        from isotope_pattern_matcher import cosine_similarity

        peaks = [(100.0, 50.0), (101.0, 20.0), (102.0, 5.0)]
        sim = cosine_similarity(peaks, peaks, mz_tolerance=0.01)
        assert abs(sim - 1.0) < 1e-9

    def test_cosine_similarity_no_overlap(self):
        from isotope_pattern_matcher import cosine_similarity

        theoretical = [(100.0, 50.0), (101.0, 20.0)]
        observed = [(200.0, 50.0), (201.0, 20.0)]
        sim = cosine_similarity(theoretical, observed, mz_tolerance=0.01)
        assert sim == 0.0

    def test_parse_peaks(self):
        from isotope_pattern_matcher import parse_peaks

        result = parse_peaks(["181.0709,100.0", "182.0742,6.7"])
        assert len(result) == 2
        assert result[0] == (181.0709, 100.0)

    def test_parse_peaks_invalid(self):
        from isotope_pattern_matcher import parse_peaks

        with pytest.raises(ValueError):
            parse_peaks(["181.0709"])
```

- [ ] **Step 5: Create README.md**

```markdown
# Isotope Pattern Generator & Matcher

Generate theoretical isotope distributions for any molecular formula and
optionally compute a cosine similarity score against observed peaks.

## Usage

```bash
python isotope_pattern_matcher.py --formula C6H12O6
python isotope_pattern_matcher.py --formula C6H12O6 --peaks 181.0709,100.0 182.0742,6.7 183.0775,0.4
```
```

- [ ] **Step 6: Run ruff and tests**

Run: `ruff check tools/metabolomics/isotope_pattern_matcher/ && PYTHONPATH=tools/metabolomics/isotope_pattern_matcher python -m pytest tools/metabolomics/isotope_pattern_matcher/tests/ -v`
Expected: Lint clean, 6 tests pass

- [ ] **Step 7: Commit**

```bash
git add tools/metabolomics/isotope_pattern_matcher/
git commit -m "Migrate isotope_pattern_matcher to per-script directory structure"
```

---

### Task 8: Migrate metabolite_feature_detection to per-script directory

**Files:**
- Create: `tools/metabolomics/metabolite_feature_detection/metabolite_feature_detection.py`
- Create: `tools/metabolomics/metabolite_feature_detection/requirements.txt`
- Create: `tools/metabolomics/metabolite_feature_detection/README.md`
- Create: `tools/metabolomics/metabolite_feature_detection/tests/conftest.py`
- Create: `tools/metabolomics/metabolite_feature_detection/tests/test_metabolite_feature_detection.py`

- [ ] **Step 1: Create directory and copy script**

```bash
mkdir -p tools/metabolomics/metabolite_feature_detection/tests
git show origin/copilot/add-agentic-scripts-for-proteomics:tools/metabolomics/metabolite_feature_detection.py > tools/metabolomics/metabolite_feature_detection/metabolite_feature_detection.py
```

- [ ] **Step 2: Create requirements.txt**

```
pyopenms
```

- [ ] **Step 3: Create conftest.py** (identical to Task 2 Step 4)

- [ ] **Step 4: Create test file with synthetic data**

Similar to Task 5, generate synthetic MSExperiment data and write to temp mzML:

```python
"""Tests for metabolite_feature_detection."""

import os
import tempfile

import pytest

from conftest import requires_pyopenms


@requires_pyopenms
class TestMetaboliteFeatureDetection:
    def test_detect_features_returns_feature_map(self):
        import pyopenms as oms
        import numpy as np
        from metabolite_feature_detection import detect_metabolite_features

        # Create a minimal synthetic experiment
        exp = oms.MSExperiment()
        for i in range(20):
            spec = oms.MSSpectrum()
            spec.setMSLevel(1)
            spec.setRT(30.0 + i * 3.0)
            mzs = np.array([180.063, 181.066, 182.070], dtype=np.float64)
            intensities = np.array([1e5, 1e4, 1e3], dtype=np.float64)
            spec.set_peaks([mzs, intensities])
            exp.addSpectrum(spec)

        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "test.mzML")
            output_path = os.path.join(tmpdir, "test.featureXML")
            oms.MzMLFile().store(input_path, exp)

            fm = detect_metabolite_features(input_path, output_path, noise_threshold=1e2)
            assert isinstance(fm, oms.FeatureMap)
            assert os.path.exists(output_path)
```

- [ ] **Step 5: Create README.md**

```markdown
# Metabolite Feature Detection

Detect small-molecule features in centroided LC-MS data using the
`FeatureFinderMetabo` pipeline. Output is written as a featureXML file.

## Usage

```bash
python metabolite_feature_detection.py --input sample.mzML
python metabolite_feature_detection.py --input sample.mzML --output features.featureXML --noise 1e5
```
```

- [ ] **Step 6: Run ruff and tests**

Run: `ruff check tools/metabolomics/metabolite_feature_detection/ && PYTHONPATH=tools/metabolomics/metabolite_feature_detection python -m pytest tools/metabolomics/metabolite_feature_detection/tests/ -v`
Expected: Lint clean, 1 test passes

- [ ] **Step 7: Commit**

```bash
git add tools/metabolomics/metabolite_feature_detection/
git commit -m "Migrate metabolite_feature_detection to per-script directory with synthetic test data"
```

---

### Task 9: Create validate-script Claude Code skill

**Files:**
- Create: `.claude/skills/validate-script.md`

- [ ] **Step 1: Create skill directory and file**

```bash
mkdir -p .claude/skills
```

Then create `.claude/skills/validate-script.md`:

```markdown
---
name: validate-script
description: Validate a pyopenms script in an isolated venv — runs ruff lint and pytest
---

# Validate Script

Validate any script in the agentomics repo by running ruff and pytest in a fresh isolated venv.

## Steps (follow exactly — rigid skill)

1. **Identify the script directory.** If the user provided a path, use it. Otherwise, ask which script to validate. The path should be `tools/<domain>/<tool_name>/`.

2. **Verify the directory structure.** Confirm it contains:
   - `<tool_name>.py`
   - `requirements.txt`
   - `tests/` directory with at least one `test_*.py` file

3. **Create a temporary venv and run validation.** Execute these commands:

   ```bash
   SCRIPT_DIR=<path-to-script-directory>
   VENV_DIR=$(mktemp -d)
   python -m venv "$VENV_DIR"
   "$VENV_DIR/bin/python" -m pip install -r "$SCRIPT_DIR/requirements.txt"
   "$VENV_DIR/bin/python" -m pip install pytest ruff
   "$VENV_DIR/bin/python" -m ruff check "$SCRIPT_DIR/"
   PYTHONPATH="$SCRIPT_DIR" "$VENV_DIR/bin/python" -m pytest "$SCRIPT_DIR/tests/" -v
   rm -rf "$VENV_DIR"
   ```

4. **Report results.** Summarize pass/fail for both ruff and pytest. If either fails, show the relevant error output so the user can fix it.

5. **Clean up.** Ensure the temporary venv is removed even if validation fails.
```

- [ ] **Step 2: Commit**

```bash
git add .claude/skills/validate-script.md
git commit -m "Add validate-script Claude Code skill for isolated venv validation"
```

---

### Task 10: Create contribute-script Claude Code skill

**Files:**
- Create: `.claude/skills/contribute-script.md`

- [ ] **Step 1: Create skill file**

```markdown
---
name: contribute-script
description: Guide creation of a new pyopenms script contribution — scaffolding through validation
---

# Contribute Script

Guide an AI agent through creating a new pyopenms CLI tool for the agentomics repo. Follow every step — this is a rigid skill.

## Prerequisites

Read `AGENTS.md` in the repo root for the full contributor guide and code patterns.

## Steps

### 1. Understand the tool

Ask the user:
- What does this tool do? What pyopenms functionality does it use?
- What gap in OpenMS/pyopenms does it fill?

### 2. Determine the domain

Ask: Is this a **proteomics** or **metabolomics** tool? If neither fits, discuss whether a new domain directory is needed.

### 3. Pick a name

Choose a descriptive snake_case name for the tool (e.g. `peptide_mass_calculator`, `isotope_pattern_matcher`). Confirm with the user.

### 4. Create a feature branch

```bash
git checkout -b add/<tool_name>
```

### 5. Scaffold the directory

```bash
mkdir -p tools/<domain>/<tool_name>/tests
```

Create these files:

**`requirements.txt`:**
```
pyopenms
```
Add any additional dependencies the script needs (one per line, no version pins).

**`tests/conftest.py`:**
```python
import sys
import os

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

try:
    import pyopenms  # noqa: F401

    HAS_PYOPENMS = True
except ImportError:
    HAS_PYOPENMS = False

requires_pyopenms = pytest.mark.skipif(not HAS_PYOPENMS, reason="pyopenms not installed")
```

### 6. Write the script

Create `tools/<domain>/<tool_name>/<tool_name>.py` following these patterns:

- Module-level docstring with description, supported features, and CLI usage examples
- pyopenms import guard:
  ```python
  try:
      import pyopenms as oms
  except ImportError:
      sys.exit("pyopenms is required. Install it with:  pip install pyopenms")
  ```
- `PROTON = 1.007276` constant where mass-to-charge calculations are needed
- Importable functions as the primary interface (with type hints and numpy-style docstrings)
- `main()` function with argparse CLI
- `if __name__ == "__main__": main()` guard

### 7. Write tests

Create `tools/<domain>/<tool_name>/tests/test_<tool_name>.py`:

- Import `requires_pyopenms` from conftest
- Decorate test classes with `@requires_pyopenms`
- Use `from <tool_name> import <function>` inside test methods
- For file-I/O scripts: generate synthetic data using pyopenms objects in test fixtures, write to `tempfile.TemporaryDirectory()`
- Cover: basic functionality, edge cases, key parameters

### 8. Write README

Create `tools/<domain>/<tool_name>/README.md` with a brief description and CLI usage examples.

### 9. Validate

Invoke the `validate-script` skill on the new script directory. Both ruff and pytest must pass.

### 10. Commit

```bash
git add tools/<domain>/<tool_name>/
git commit -m "Add <tool_name>: <brief description>"
```
```

- [ ] **Step 2: Commit**

```bash
git add .claude/skills/contribute-script.md
git commit -m "Add contribute-script Claude Code skill for guided new tool creation"
```

---

### Task 11: Create AGENTS.md

**Files:**
- Create: `AGENTS.md`

- [ ] **Step 1: Create AGENTS.md**

```markdown
# AGENTS.md — AI Contributor Guide

This file instructs AI agents (Claude Code, GitHub Copilot, Cursor, Gemini, etc.) how to contribute scripts to the agentomics repository.

## Project Purpose

Agentomics is a collection of standalone CLI tools built with [pyopenms](https://pyopenms.readthedocs.io/) for proteomics and metabolomics workflows. These tools fill gaps not yet covered by OpenMS/pyopenms. All code in this repo is written by AI agents.

## Contribution Requirements

Every script must be a **self-contained directory** under `tools/<domain>/<tool_name>/`:

```
tools/<domain>/<tool_name>/
├── <tool_name>.py        # The tool itself
├── requirements.txt      # pyopenms + any script-specific deps (no version pins)
├── README.md             # Brief description + CLI usage examples
└── tests/
    ├── conftest.py       # Shared test config (see below)
    └── test_<tool_name>.py
```

### Rules

- `<domain>` is `proteomics` or `metabolomics`
- `requirements.txt` always includes `pyopenms` with no version pin — builds against latest
- No cross-script imports — each script is fully independent
- No `__init__.py` files — these are NOT Python packages
- No scripts that duplicate functionality already in OpenMS/pyopenms

## Code Patterns

### Script structure

Every script must have:

1. **Module docstring** with description, features, and usage examples
2. **pyopenms import guard:**
   ```python
   import sys
   try:
       import pyopenms as oms
   except ImportError:
       sys.exit("pyopenms is required. Install it with:  pip install pyopenms")
   ```
3. **Importable functions** as the primary interface (with type hints and numpy-style docstrings)
4. **`main()` function** with argparse CLI
5. **`if __name__ == "__main__": main()`** guard
6. **`PROTON = 1.007276`** constant where mass-to-charge calculations are needed

### Test structure

Every `tests/conftest.py` must contain:

```python
import sys
import os

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

try:
    import pyopenms  # noqa: F401
    HAS_PYOPENMS = True
except ImportError:
    HAS_PYOPENMS = False

requires_pyopenms = pytest.mark.skipif(not HAS_PYOPENMS, reason="pyopenms not installed")
```

Test files:
- Decorate test classes with `@requires_pyopenms` from conftest
- Import script functions inside test methods: `from <tool_name> import <function>`
- For file-I/O scripts: generate synthetic data using pyopenms objects, write to `tempfile.TemporaryDirectory()`

## Validation

Every script must pass validation in an **isolated venv** before it can be merged. Run these commands from the repo root:

```bash
SCRIPT_DIR=tools/<domain>/<tool_name>
VENV_DIR=$(mktemp -d)
python -m venv "$VENV_DIR"
"$VENV_DIR/bin/python" -m pip install -r "$SCRIPT_DIR/requirements.txt"
"$VENV_DIR/bin/python" -m pip install pytest ruff
"$VENV_DIR/bin/python" -m ruff check "$SCRIPT_DIR/"
PYTHONPATH="$SCRIPT_DIR" "$VENV_DIR/bin/python" -m pytest "$SCRIPT_DIR/tests/" -v
rm -rf "$VENV_DIR"
```

Both ruff and pytest must pass with zero errors.

## Linting

Ruff is configured in `ruff.toml` at the repo root:
- Line length: 120
- Rules: E (pycodestyle errors), F (pyflakes), W (pycodestyle warnings), I (isort)

## What NOT to Do

- Do not add cross-script imports
- Do not add dependencies to a shared/root requirements file
- Do not create scripts that duplicate existing pyopenms CLI tools or OpenMS TOPP tools
- Do not pin pyopenms to a specific version
- Do not add `__init__.py` files
```

- [ ] **Step 2: Commit**

```bash
git add AGENTS.md
git commit -m "Add AGENTS.md platform-agnostic AI contributor guide"
```

---

### Task 12: Create GitHub Actions CI workflow

**Files:**
- Create: `.github/workflows/validate.yml`

- [ ] **Step 1: Create workflow file**

```yaml
name: Validate Scripts

on:
  pull_request:
    paths:
      - 'tools/**'

jobs:
  detect-changes:
    runs-on: ubuntu-latest
    outputs:
      matrix: ${{ steps.detect.outputs.matrix }}
      has_changes: ${{ steps.detect.outputs.has_changes }}
    steps:
      - uses: actions/checkout@v4
        with:
          fetch-depth: 0

      - id: detect
        name: Detect changed script directories
        run: |
          # Note: github.base_ref is only available on pull_request events
          # Find all script directories that changed in this PR
          CHANGED=$(git diff --name-only origin/${{ github.base_ref }}...HEAD -- 'tools/' \
            | grep -oP 'tools/[^/]+/[^/]+/' \
            | sort -u \
            | jq -R -s -c 'split("\n") | map(select(length > 0))')

          if [ "$CHANGED" = "[]" ] || [ -z "$CHANGED" ]; then
            echo "has_changes=false" >> "$GITHUB_OUTPUT"
            echo "matrix=[]" >> "$GITHUB_OUTPUT"
          else
            echo "has_changes=true" >> "$GITHUB_OUTPUT"
            echo "matrix=$CHANGED" >> "$GITHUB_OUTPUT"
          fi

  validate:
    needs: detect-changes
    if: needs.detect-changes.outputs.has_changes == 'true'
    runs-on: ubuntu-latest
    strategy:
      fail-fast: false
      matrix:
        script_dir: ${{ fromJson(needs.detect-changes.outputs.matrix) }}
    name: Validate ${{ matrix.script_dir }}
    steps:
      - uses: actions/checkout@v4

      - uses: actions/setup-python@v5
        with:
          python-version: '3.11'

      - name: Create venv and install dependencies
        run: |
          python -m venv /tmp/validate_venv
          /tmp/validate_venv/bin/python -m pip install -r ${{ matrix.script_dir }}requirements.txt
          /tmp/validate_venv/bin/python -m pip install pytest ruff

      - name: Lint with ruff
        run: |
          /tmp/validate_venv/bin/python -m ruff check ${{ matrix.script_dir }}

      - name: Run tests
        run: |
          PYTHONPATH=${{ matrix.script_dir }} /tmp/validate_venv/bin/python -m pytest ${{ matrix.script_dir }}tests/ -v
```

- [ ] **Step 2: Commit**

```bash
mkdir -p .github/workflows
git add .github/workflows/validate.yml
git commit -m "Add GitHub Actions CI workflow for per-script isolated validation"
```

---

### Task 13: Update CLAUDE.md and README.md

**Files:**
- Modify: `CLAUDE.md`
- Modify: `README.md`

- [ ] **Step 1: Update CLAUDE.md**

Replace the current content with updated version reflecting the new structure:

```markdown
# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Purpose

Agentomics is a collection of standalone CLI tools built with [pyopenms](https://pyopenms.readthedocs.io/) for proteomics and metabolomics workflows. These tools fill gaps not yet covered by OpenMS/pyopenms. All code in this repo is agentic-only development — written entirely by AI agents.

## Commands

```bash
# Install dependencies for a specific script
pip install -r tools/proteomics/peptide_mass_calculator/requirements.txt

# Lint a specific script
ruff check tools/proteomics/peptide_mass_calculator/

# Run tests for a specific script
PYTHONPATH=tools/proteomics/peptide_mass_calculator python -m pytest tools/proteomics/peptide_mass_calculator/tests/ -v

# Lint all scripts
ruff check tools/

# Run all tests (requires pyopenms installed)
find scripts -name 'tests' -type d -exec sh -c 'PYTHONPATH=$(dirname {}) python -m pytest {} -v' \;
```

## Architecture

### Per-Script Directory Structure

Each script is a self-contained directory under `tools/<domain>/<tool_name>/`:

```
tools/<domain>/<tool_name>/
├── <tool_name>.py        # The tool (importable functions + argparse CLI)
├── requirements.txt      # pyopenms + script-specific deps
├── README.md             # Usage examples
└── tests/
    ├── conftest.py       # requires_pyopenms marker + sys.path setup
    └── test_<tool_name>.py
```

Domains: `proteomics/`, `metabolomics/`

### Key Patterns

- pyopenms import wrapped in try/except with user-friendly error message
- Mass-to-charge: `(mass + charge * PROTON) / charge` with `PROTON = 1.007276`
- Every script has dual interface: importable functions + argparse CLI + `__main__` guard
- Tests use `@requires_pyopenms` skip marker from conftest.py
- File-I/O scripts use synthetic test data generated with pyopenms objects

## Contributing

See `AGENTS.md` for the full AI contributor guide. Two Claude Code skills are available:

- **`contribute-script`** — guided workflow for adding a new script
- **`validate-script`** — validate any script in an isolated venv (ruff + pytest)
```

- [ ] **Step 2: Update README.md**

```markdown
# agentomics

A repository of agentic-created tools using [pyopenms](https://pyopenms.readthedocs.io/) for proteomics and metabolomics.

All code in this repo is written by AI agents. See [AGENTS.md](AGENTS.md) for the contributor guide.

## Requirements

```bash
pip install pyopenms
```

## Scripts

### Proteomics

| Script | Description |
|--------|-------------|
| [`peptide_mass_calculator`](tools/proteomics/peptide_mass_calculator/) | Monoisotopic/average masses and b/y fragment ions for peptide sequences |
| [`protein_digest`](tools/proteomics/protein_digest/) | In-silico enzymatic protein digestion |
| [`spectrum_file_info`](tools/proteomics/spectrum_file_info/) | Summary statistics for mzML files |
| [`feature_detection_proteomics`](tools/proteomics/feature_detection_proteomics/) | Peptide feature detection from LC-MS/MS data |

### Metabolomics

| Script | Description |
|--------|-------------|
| [`mass_accuracy_calculator`](tools/metabolomics/mass_accuracy_calculator/) | m/z mass accuracy (ppm error) for sequences or formulas |
| [`isotope_pattern_matcher`](tools/metabolomics/isotope_pattern_matcher/) | Theoretical isotope distributions and cosine similarity scoring |
| [`metabolite_feature_detection`](tools/metabolomics/metabolite_feature_detection/) | Metabolite feature detection from LC-MS data |

## Validation

Each script is validated in an isolated venv. See [AGENTS.md](AGENTS.md) for validation commands.

## License

BSD 3-Clause — see [LICENSE](LICENSE).
```

- [ ] **Step 3: Commit**

```bash
git add CLAUDE.md README.md
git commit -m "Update CLAUDE.md and README.md for per-script directory structure"
```

---

### Task 14: Final validation — run all scripts through isolated validation

- [ ] **Step 1: Validate all 7 scripts**

Run the validation pipeline on each script directory. For each, execute:

```bash
SCRIPT_DIR=tools/<domain>/<tool_name>
VENV_DIR=$(mktemp -d)
python -m venv "$VENV_DIR"
"$VENV_DIR/bin/python" -m pip install -r "$SCRIPT_DIR/requirements.txt"
"$VENV_DIR/bin/python" -m pip install pytest ruff
"$VENV_DIR/bin/python" -m ruff check "$SCRIPT_DIR/"
PYTHONPATH="$SCRIPT_DIR" "$VENV_DIR/bin/python" -m pytest "$SCRIPT_DIR/tests/" -v
rm -rf "$VENV_DIR"
```

Run for each:
1. `tools/proteomics/peptide_mass_calculator`
2. `tools/proteomics/protein_digest`
3. `tools/proteomics/spectrum_file_info`
4. `tools/proteomics/feature_detection_proteomics`
5. `tools/metabolomics/mass_accuracy_calculator`
6. `tools/metabolomics/isotope_pattern_matcher`
7. `tools/metabolomics/metabolite_feature_detection`

Expected: All 7 pass ruff lint and all tests pass (or skip with `pyopenms not installed`).

- [ ] **Step 2: Fix any failures and recommit**

If any script fails lint or tests, fix the issue and create a new commit with the fix.

- [ ] **Step 3: Final commit if any fixes were needed**

```bash
git add -u
git commit -m "Fix validation issues found during final check"
```
