# pyopenms Gap Coverage Tools — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement 50 new pyopenms CLI tools covering signal processing, alignment, quantification, inference, DIA, deconvolution, search, and metabolomics preprocessing.

**Architecture:** Each tool is a self-contained directory under `tools/<domain>/<topic>/<tool_name>/` with a click CLI, importable functions, and pytest tests. All tools follow the identical scaffolding pattern. 12 independent topic areas can be fully parallelized.

**Tech Stack:** Python 3.11+, pyopenms, click, pytest, ruff

**Spec:** `docs/superpowers/specs/2026-03-25-pyopenms-gap-tools-design.md`

---

## Scaffolding Template (referenced by all tasks)

Every tool uses this exact structure. Replace `{DOMAIN}`, `{TOPIC}`, `{TOOL}`, etc.

**Directory:** `tools/{DOMAIN}/{TOPIC}/{TOOL}/`

**`requirements.txt`:**
```
pyopenms
click
```

**`tests/conftest.py`:**
```python
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
```

**`{TOOL}.py` skeleton:**
```python
"""
{TITLE}
{'=' * len(TITLE)}
{DESCRIPTION}

Usage
-----
    python {TOOL}.py {CLI_EXAMPLE}
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def {FUNCTION}({PARAMS}) -> {RETURN_TYPE}:
    """{DOCSTRING}"""
    # implementation


@click.command()
@click.option(...)
def main(...):
    """{CLI_HELP}"""
    result = {FUNCTION}(...)
    print(f"...")


if __name__ == "__main__":
    main()
```

**`tests/test_{TOOL}.py` skeleton:**
```python
"""Tests for {TOOL}."""

import pytest

pytest.importorskip("pyopenms")


class Test{ClassName}:
    def test_{BASIC_CASE}(self):
        from {TOOL} import {FUNCTION}
        # test body
```

**`README.md` skeleton:**
```markdown
# {TITLE}

{DESCRIPTION}

## Installation

\`\`\`bash
pip install -r requirements.txt
\`\`\`

## Usage

\`\`\`bash
python {TOOL}.py {CLI_EXAMPLE}
\`\`\`
```

**Validation per tool:**
```bash
ruff check tools/{DOMAIN}/{TOPIC}/{TOOL}/
PYTHONPATH=tools/{DOMAIN}/{TOPIC}/{TOOL} python -m pytest tools/{DOMAIN}/{TOPIC}/{TOOL}/tests/ -v
```

---

## Task 0: Update topic directories in AGENTS.md and CLAUDE.md

**Files:**
- Modify: `AGENTS.md` (proteomics topics list, ~line 25)
- Modify: `CLAUDE.md` (proteomics and metabolomics topics lists, ~lines 50-52)

- [ ] **Step 1: Add new proteomics topics to AGENTS.md**

Add to the proteomics topics line: `signal_processing/`, `map_alignment/`, `quantification/`, `inference/`, `dia_analysis/`

- [ ] **Step 2: Add new metabolomics topics to AGENTS.md**

Add to the metabolomics topics line: `identification/`, `preprocessing/`

- [ ] **Step 3: Add same topics to CLAUDE.md**

Mirror the changes in CLAUDE.md proteomics and metabolomics topic lists.

- [ ] **Step 4: Commit**

```bash
git add AGENTS.md CLAUDE.md
git commit -m "Add new topic directories for pyopenms gap tools"
```

---

## Task 1: Signal Processing — Smoothing & Baseline (4 tools)

**Topic:** `proteomics/signal_processing/`

**Tools to implement:**

### 1.1 `gaussian_smoother`
- **Wraps:** `oms.GaussFilter()`
- **Key code:** `gf = oms.GaussFilter(); params = gf.getDefaults(); params.setValue("gaussian_width", gaussian_width); gf.setParameters(params); gf.filterExperiment(exp)`
- **CLI:** `--input run.mzML --output smoothed.mzML --gaussian-width 0.2`
- **Function:** `smooth_experiment(input_path, output_path, gaussian_width=0.2)` → int
- **Tests:** Create synthetic 5-scan mzML with noisy Gaussian profiles via `spec.set_peaks()`, smooth, verify intensity variance reduced

### 1.2 `savitzky_golay_smoother`
- **Wraps:** `oms.SavitzkyGolayFilter()`
- **Key code:** `sgf = oms.SavitzkyGolayFilter(); params = sgf.getDefaults(); params.setValue("frame_length", frame_length); params.setValue("polynomial_order", polynomial_order); sgf.setParameters(params); sgf.filterExperiment(exp)`
- **CLI:** `--input run.mzML --output smoothed.mzML --frame-length 11 --polynomial-order 3`
- **Function:** `smooth_experiment(input_path, output_path, frame_length=11, polynomial_order=3)` → int
- **Tests:** Synthetic profile data, verify output spectrum count matches input, peaks preserved

### 1.3 `baseline_corrector`
- **Wraps:** `oms.MorphologicalFilter()`
- **Key code:** `mf = oms.MorphologicalFilter(); params = mf.getDefaults(); params.setValue("struc_elem_length", struct_element_length); mf.setParameters(params); mf.filterExperiment(exp)`
- **CLI:** `--input run.mzML --output corrected.mzML --struct-element-length 3.0`
- **Function:** `correct_baseline(input_path, output_path, struct_element_length=3.0)` → int
- **Tests:** Synthetic spectrum with uniform baseline offset + peaks, verify baseline reduced while peaks persist

### 1.4 `spectrum_normalizer`
- **Wraps:** `oms.Normalizer()`
- **Key code:** `norm = oms.Normalizer(); params = norm.getDefaults(); params.setValue("method", method); norm.setParameters(params); norm.filterPeakMap(exp)`
- **CLI:** `--input run.mzML --output normalized.mzML --method to_one`
- **Function:** `normalize_experiment(input_path, output_path, method="to_one")` → int
- **Tests:** Synthetic spectrum, verify max intensity = 1.0 (`to_one`) or sum ≈ 1.0 (`to_TIC`)

**Steps per tool (repeat for each):**

- [ ] **Step 1:** Scaffold directory: `mkdir -p tools/proteomics/signal_processing/{TOOL}/tests`
- [ ] **Step 2:** Write `requirements.txt`, `tests/conftest.py`, `README.md`
- [ ] **Step 3:** Write test file with 2-4 test cases using synthetic data
- [ ] **Step 4:** Run tests to verify they fail: `PYTHONPATH=tools/proteomics/signal_processing/{TOOL} python -m pytest tools/proteomics/signal_processing/{TOOL}/tests/ -v`
- [ ] **Step 5:** Implement the tool
- [ ] **Step 6:** Run tests to verify they pass
- [ ] **Step 7:** Run `ruff check tools/proteomics/signal_processing/{TOOL}/`
- [ ] **Step 8:** Commit: `git add tools/proteomics/signal_processing/{TOOL}/ && git commit -m "Add {TOOL}: {description}"`

---

## Task 2: Signal Processing — Noise Filters (3 tools)

**Topic:** `proteomics/signal_processing/`

### 2.1 `noise_filter_window`
- **Wraps:** `oms.WindowMower()`
- **Key code:** `wm = oms.WindowMower(); params = wm.getDefaults(); params.setValue("windowsize", window_size); params.setValue("peakcount", peak_count); wm.setParameters(params); wm.filterPeakMap(exp)`
- **CLI:** `--input run.mzML --output filtered.mzML --window-size 50.0 --peak-count 3`
- **Function:** `filter_window(input_path, output_path, window_size=50.0, peak_count=3)` → int
- **Tests:** Synthetic spectrum with 20 peaks, filter with peak_count=5, verify ≤5 peaks per window

### 2.2 `noise_filter_threshold`
- **Wraps:** `oms.ThresholdMower()`
- **Key code:** `tm = oms.ThresholdMower(); params = tm.getDefaults(); params.setValue("threshold", threshold); tm.setParameters(params); tm.filterPeakMap(exp)`
- **CLI:** `--input run.mzML --output filtered.mzML --threshold 100.0`
- **Function:** `filter_threshold(input_path, output_path, threshold=100.0)` → int
- **Tests:** Synthetic spectrum with peaks at intensities 50, 150, 500; threshold=100; verify only 150 and 500 remain

### 2.3 `noise_filter_nlargest`
- **Wraps:** `oms.NLargest()`
- **Key code:** `nl = oms.NLargest(); params = nl.getDefaults(); params.setValue("n", n); nl.setParameters(params); nl.filterPeakMap(exp)`
- **CLI:** `--input run.mzML --output filtered.mzML --n 100`
- **Function:** `filter_nlargest(input_path, output_path, n=100)` → int
- **Tests:** Spectrum with 200 peaks, n=50, verify output has ≤50 peaks per spectrum

**Steps:** Same as Task 1 (scaffold → test → implement → lint → commit per tool)

---

## Task 3: Signal Processing — S/N, Peak Picking, Resampling (5 tools)

**Topic:** `proteomics/signal_processing/`

### 3.1 `signal_to_noise_estimator`
- **Wraps:** `oms.SignalToNoiseEstimatorMedian()`
- **Key code:** `sne = oms.SignalToNoiseEstimatorMedian(); sne.init(spectrum); sn = sne.getSignalToNoise(i)` for each peak index
- **CLI:** `--input run.mzML --output sn_report.tsv`
- **Function:** `estimate_sn(input_path, output_path)` → int (spectra analyzed)
- **Output TSV columns:** spectrum_index, mz, intensity, sn_ratio
- **Tests:** Synthetic spectrum with high-intensity peak among low noise, verify S/N > 10 for strong peak

### 3.2 `peak_picker_hires`
- **Wraps:** `oms.PeakPickerHiRes()`
- **Key code:** `pp = oms.PeakPickerHiRes(); params = pp.getDefaults(); params.setValue("signal_to_noise", sn); pp.setParameters(params); pp.pickExperiment(exp_in, exp_out)`
- **CLI:** `--input profile.mzML --output centroid.mzML --signal-to-noise 1.0`
- **Function:** `pick_peaks(input_path, output_path, signal_to_noise=1.0)` → int
- **Tests:** Synthetic profile mzML (dense m/z points forming Gaussian peaks), verify output has fewer points, peaks at correct m/z

### 3.3 `chromatogram_peak_picker`
- **Wraps:** `oms.PeakPickerChromatogram()`
- **Key code:** `pp = oms.PeakPickerChromatogram(); pp.pickChromatogram(chrom_in, chrom_out)`
- **CLI:** `--input run.mzML --output picked.mzML`
- **Function:** `pick_chromatogram_peaks(input_path, output_path)` → int
- **Tests:** Synthetic mzML with chromatograms having Gaussian RT profiles, verify peak detection

### 3.4 `peak_type_estimator`
- **Wraps:** `oms.PeakTypeEstimator()`
- **Key code:** `pte = oms.PeakTypeEstimator(); peak_type = pte.estimateType(spectrum)`
- **CLI:** `--input run.mzML --output report.tsv`
- **Function:** `estimate_peak_type(input_path, output_path)` → dict with type and counts
- **Output TSV columns:** spectrum_index, ms_level, peak_type
- **Tests:** One synthetic profile spectrum (dense points), one centroided (sparse), verify correct classification

### 3.5 `linear_resampler`
- **Wraps:** `oms.LinearResampler()`
- **Key code:** `lr = oms.LinearResampler(); params = lr.getDefaults(); params.setValue("spacing", spacing); lr.setParameters(params); lr.rasterExperiment(exp)`
- **CLI:** `--input run.mzML --output resampled.mzML --spacing 0.01`
- **Function:** `resample_experiment(input_path, output_path, spacing=0.01)` → int
- **Tests:** Synthetic spectrum, verify output m/z values are uniformly spaced

**Steps:** Same as Task 1

---

## Task 4: Map Alignment (3 tools)

**Topic:** `proteomics/map_alignment/`

### 4.1 `rt_alignment_pose_clustering`
- **Wraps:** `oms.MapAlignmentAlgorithmPoseClustering()`
- **Key code:**
```python
aligner = oms.MapAlignmentAlgorithmPoseClustering()
aligner.setReference(reference_map)
trafo = oms.TransformationDescription()
aligner.align(input_map, trafo)
oms.MapAlignmentTransformer.transformRetentionTimes(input_map, trafo)
```
- **CLI:** `--reference ref.featureXML --input sample.featureXML --output aligned.featureXML`
- **Function:** `align_pose_clustering(reference_path, input_path, output_path)` → dict
- **Tests:** Create 2 synthetic FeatureMaps (one with RT+10.0 offset), align, verify RT difference < 1.0

### 4.2 `rt_alignment_identification`
- **Wraps:** `oms.MapAlignmentAlgorithmIdentification()`
- **Key code:** Similar pattern; uses idXML (ProteinIdentification + PeptideIdentification lists)
- **CLI:** `--reference ref.idXML --input sample.idXML --output aligned.idXML --model linear`
- **Function:** `align_identification(reference_path, input_path, output_path, model="linear")` → dict
- **Tests:** Synthetic idXML with peptide hits at shifted RTs, verify alignment

### 4.3 `rt_alignment_kd`
- **Wraps:** `oms.MapAlignmentAlgorithmKD()`
- **CLI:** `--inputs run1.featureXML run2.featureXML --output-dir aligned/ --model lowess`
- **Function:** `align_kd(input_paths, output_dir, model="lowess")` → dict
- **Tests:** Two synthetic FeatureMaps with systematic RT offset

**Steps:** Same as Task 1

---

## Task 5: Feature Linking (3 tools)

**Topic:** `proteomics/map_alignment/`

### 5.1 `feature_linker_qt`
- **Wraps:** `oms.FeatureGroupingAlgorithmQT()`
- **Key code:**
```python
linker = oms.FeatureGroupingAlgorithmQT()
params = linker.getDefaults()
params.setValue("distance_RT:max_difference", rt_tol)
params.setValue("distance_MZ:max_difference", mz_tol)
params.setValue("distance_MZ:unit", "Da")
linker.setParameters(params)
consensus = oms.ConsensusMap()
file_descriptions = consensus.getColumnHeaders()
# set up file descriptions for each input
linker.group(feature_maps, consensus)
```
- **CLI:** `--inputs run1.featureXML run2.featureXML --output consensus.consensusXML --rt-tol 30.0 --mz-tol 0.01`
- **Function:** `link_features_qt(input_paths, output_path, rt_tol=30.0, mz_tol=0.01)` → int
- **Tests:** 2 synthetic FeatureMaps with 3 features each (same m/z, close RT), verify 3 consensus features

### 5.2 `feature_linker_kd`
- **Wraps:** `oms.FeatureGroupingAlgorithmKD()`
- **Key code:** Same pattern as QT with `FeatureGroupingAlgorithmKD`
- **CLI:** `--inputs run1.featureXML run2.featureXML --output consensus.consensusXML --rt-tol 30.0 --mz-tol 0.01`
- **Function:** `link_features_kd(input_paths, output_path, rt_tol=30.0, mz_tol=0.01)` → int
- **Tests:** Same as QT

### 5.3 `feature_linker_labeled`
- **Wraps:** `oms.FeatureGroupingAlgorithmLabeled()`
- **CLI:** `--input run.featureXML --output consensus.consensusXML --rt-estimate 1000.0 --mz-estimate 5.0`
- **Function:** `link_labeled_features(input_path, output_path, rt_estimate=1000.0, mz_estimate=5.0)` → int
- **Tests:** Synthetic FeatureMap with light/heavy pairs (Δm ≈ 8 Da for SILAC K+8), verify pairing

**Steps:** Same as Task 1

---

## Task 6: Quantification — Isobaric (3 tools)

**Topic:** `proteomics/quantification/`

### 6.1 `isobaric_channel_extractor`
- **Wraps:** `oms.IsobaricChannelExtractor()`
- **Key code:**
```python
method_map = {"tmt6plex": oms.TMTSixPlexQuantitationMethod(), "tmt10plex": oms.TMTTenPlexQuantitationMethod(), ...}
quant_method = method_map[method]
extractor = oms.IsobaricChannelExtractor(quant_method)
consensus = oms.ConsensusMap()
extractor.extractChannels(exp, consensus)
```
- **CLI:** `--input run.mzML --method tmt6plex --output quant.consensusXML`
- **Function:** `extract_channels(input_path, method, output_path)` → int
- **Tests:** Synthetic MS2 with peaks at TMT 6-plex reporter m/z values, verify extraction

### 6.2 `isobaric_isotope_corrector`
- **Wraps:** `oms.IsobaricIsotopeCorrector()`
- **Key code:** `corrector = oms.IsobaricIsotopeCorrector(); corrector.correctIsotopicImpurities(consensus, quant_method)`
- **CLI:** `--input quant.consensusXML --method tmt6plex --output corrected.consensusXML`
- **Function:** `correct_isotope_impurities(input_path, method, output_path)` → int
- **Tests:** Consensus map with known impurities, verify correction changes intensities

### 6.3 `isobaric_quantifier`
- **Wraps:** `oms.IsobaricQuantifier()`
- **Key code:** `quant = oms.IsobaricQuantifier(quant_method); consensus_out = oms.ConsensusMap(); quant.quantify(consensus_in, consensus_out)`
- **CLI:** `--input quant.consensusXML --method tmt6plex --output quantified.consensusXML`
- **Function:** `quantify_isobaric(input_path, method, output_path)` → int
- **Tests:** Synthetic consensus map, verify quantification produces non-zero values

**Steps:** Same as Task 1

---

## Task 7: Quantification — Multiplex, Normalization, LFQ (5 tools)

**Topic:** `proteomics/quantification/`

### 7.1 `multiplex_feature_finder`
- **Wraps:** `oms.FeatureFinderMultiplexAlgorithm()`
- **CLI:** `--input run.mzML --output features.featureXML --labels light:0 heavy:8.0142`
- **Function:** `find_multiplex_features(input_path, output_path, labels)` → int
- **Tests:** Synthetic mzML with light/heavy peptide pairs separated by 8.0142 Da

### 7.2 `consensus_map_normalizer`
- **Wraps:** `oms.ConsensusMapNormalizerAlgorithmMedian`, `oms.ConsensusMapNormalizerAlgorithmQuantile`
- **Key code:** `oms.ConsensusMapNormalizerAlgorithmMedian.normalizeMaps(consensus, "median", "", [])` — note: check `normalizeMaps` signature at runtime for correct enum/int normalization method parameter
- **CLI:** `--input consensus.consensusXML --output normalized.consensusXML --method median`
- **Function:** `normalize_consensus(input_path, output_path, method="median")` → int
- **Tests:** Synthetic consensus map with 2x intensity bias in one channel, verify normalization

### 7.3 `id_mapper`
- **Wraps:** `oms.IDMapper()`
- **Key code:** `mapper = oms.IDMapper(); params = mapper.getDefaults(); params.setValue("rt_tolerance", rt_tol); params.setValue("mz_tolerance", mz_tol); mapper.setParameters(params); mapper.annotate(feature_map, peptide_ids, protein_ids)`
- **CLI:** `--features features.featureXML --ids peptides.idXML --output annotated.featureXML --rt-tol 5.0 --mz-tol 10.0`
- **Function:** `map_ids_to_features(features_path, ids_path, output_path, rt_tol=5.0, mz_tol=10.0)` → int
- **Tests:** Synthetic feature at RT=100, m/z=500 + peptide hit at RT=100.5, m/z=500.001, verify annotation

### 7.4 `peptide_protein_quantifier`
- **Wraps:** `oms.PeptideAndProteinQuant()`
- **CLI:** `--input annotated.consensusXML --output protein_quant.csv`
- **Function:** `quantify_proteins(input_path, output_path)` → int
- **Tests:** Synthetic consensus map with peptide annotations from 2 proteins, verify rollup

### 7.5 `id_conflict_resolver`
- **Wraps:** `oms.IDConflictResolverAlgorithm()`
- **Key code:** `resolver = oms.IDConflictResolverAlgorithm(); resolver.resolve(feature_map)`
- **CLI:** `--input features.featureXML --output resolved.featureXML`
- **Function:** `resolve_id_conflicts(input_path, output_path)` → int
- **Tests:** Feature with 2 peptide hits (scores 0.9 and 0.5), verify only best retained

**Steps:** Same as Task 1

---

## Task 8: Inference / FDR (7 tools)

**Topic:** `proteomics/inference/`

### 8.1 `fdr_estimator`
- **Wraps:** `oms.FalseDiscoveryRate()`
- **Key code:** `fdr = oms.FalseDiscoveryRate(); fdr.apply(peptide_ids)` (also `fdr.apply(protein_ids)` for protein-level FDR)
- **Function:** `estimate_fdr(input_path, output_path, protein_level=False)` → dict
- **Tests:** Synthetic idXML with 5 target (high score) + 5 decoy (low score) PSMs, verify q-values assigned

### 8.2 `posterior_error_probability`
- **Wraps:** `oms.PosteriorErrorProbabilityModel()`
- **Function:** `estimate_pep(input_path, output_path)` → int
- **Tests:** Synthetic PSMs with bimodal score distribution, verify PEP values between 0-1

### 8.3 `protein_inference_basic`
- **Wraps:** `oms.BasicProteinInferenceAlgorithm()`
- **Key code:** `algo = oms.BasicProteinInferenceAlgorithm(); algo.run(peptide_ids, protein_ids)`
- **Function:** `infer_proteins(input_path, output_path)` → int
- **Tests:** 3 peptides mapping to 2 proteins (one shared), verify both proteins inferred

### 8.4 `protein_inference_bayesian`
- **Wraps:** `oms.BayesianProteinInferenceAlgorithm()`
- **Function:** `infer_proteins_bayesian(input_path, output_path)` → int
- **Tests:** Same as basic but verify posterior probabilities

### 8.5 `peptide_indexer`
- **Wraps:** `oms.PeptideIndexing()`
- **Key code:** `indexer = oms.PeptideIndexing(); indexer.run(fasta_db, protein_ids, peptide_ids)`
- **Function:** `index_peptides(ids_path, fasta_path, output_path)` → int
- **Tests:** FASTA with known protein, idXML with peptide from that protein, verify protein reference added

### 8.6 `id_filter`
- **Wraps:** `oms.IDFilter()`
- **Key code:** `oms.IDFilter.filterHitsByScore(peptide_ids, threshold)`, `oms.IDFilter.removeDecoyHits(peptide_ids)`
- **Function:** `filter_ids(input_path, output_path, score_threshold=0.05, remove_decoys=False)` → dict
- **Tests:** 10 PSMs with varied scores + decoys, filter at threshold, verify correct count retained

### 8.7 `consensus_id`
- **Wraps:** `oms.ConsensusIDAlgorithmAverage`, `oms.ConsensusIDAlgorithmBest`, `oms.ConsensusIDAlgorithmRanks`
- **Function:** `consensus_id(input_paths, output_path, algorithm="best")` → int
- **Tests:** Two idXML with overlapping peptide hits, verify consensus merging

**Steps:** Same as Task 1

---

## Task 9: Deconvolution (3 tools)

### 9.1 `deisotoper` — `proteomics/spectrum_analysis/`
- **Wraps:** `oms.Deisotoper()`
- **Key code:** `oms.Deisotoper.deisotopeAndSingleCharge(spectrum, fragment_tolerance, False, min_charge, max_charge, False, 3, 10, False, True)`
- **Function:** `deisotope(input_path, output_path, fragment_tolerance=0.1, min_charge=1, max_charge=5)` → int
- **Tests:** Synthetic spectrum with isotope envelope (m/z, m/z+1.003, m/z+2.006), verify reduced to monoisotopic

### 9.2 `flash_deconv` — `proteomics/spectrum_analysis/`
- **Wraps:** `oms.FLASHDeconvAlgorithm()`
- **Note:** Complex implementation — many parameters, non-trivial test data generation
- **Function:** `deconvolve_intact(input_path, output_path, min_mass=5000, max_mass=100000)` → int
- **Tests:** Synthetic multiply-charged peaks at m/z = (M + z*1.007)/z for M=10000, z=5..20, verify M recovered within tolerance

### 9.3 `metabolite_feature_deconvolution` — `metabolomics/feature_processing/`
- **Wraps:** `oms.MetaboliteFeatureDeconvolution()`
- **Key code:** `mfd = oms.MetaboliteFeatureDeconvolution(); out_fm = oms.FeatureMap(); out_cm = oms.ConsensusMap(); mfd.compute(in_fm, out_fm, out_cm)`
- **Function:** `deconvolve_adducts(input_path, output_path)` → int
- **Tests:** Synthetic features at [M+H]+ (m/z=181.07), [M+Na]+ (m/z=203.05), [M+K]+ (m/z=219.03) for glucose, verify grouping

**Steps:** Same as Task 1

---

## Task 10: OpenSWATH / DIA (6 tools)

**Topic:** `proteomics/dia_analysis/`

### 10.1 `chromatogram_extractor`
- **Wraps:** `oms.ChromatogramExtractor()`
- **Function:** `extract_chromatograms(input_path, transitions_path, output_path, mz_tol=10.0)` → int
- **Tests:** Synthetic DIA mzML with known m/z peaks, TSV transitions file, verify XIC extraction

### 10.2 `mrm_feature_finder`
- **Wraps:** `oms.MRMFeatureFinderScoring()`
- **Function:** `find_mrm_features(input_path, transitions_path, output_path)` → int
- **Tests:** Synthetic chromatograms with Gaussian peaks, verify feature detection

### 10.3 `dia_scorer`
- **Wraps:** `oms.DIAScoring()`
- **Function:** `score_dia(input_path, transitions_path, output_path)` → int
- **Tests:** Synthetic DIA spectrum with transitions, verify score output

### 10.4 `mrm_transition_group_picker`
- **Wraps:** `oms.MRMTransitionGroupPicker()`
- **Function:** `pick_transition_groups(input_path, output_path)` → int
- **Tests:** Synthetic transition group chromatograms, verify picking

### 10.5 `mrm_rt_normalizer`
- **Wraps:** `oms.MRMRTNormalizer()` (static utility methods: `removeOutliersIterative`, `removeOutliersRANSAC`, `computeBinnedCoverage`)
- **Note:** `MRMRTNormalizer` provides only static outlier-removal helpers, not a complete normalization workflow. The tool must: (1) read feature pairs (measured RT vs reference iRT), (2) apply outlier removal via `removeOutliersRANSAC` or `removeOutliersIterative`, (3) fit a `TransformationModelLinear` or similar, (4) apply the transformation. This is more complex than other tools.
- **Function:** `normalize_mrm_rt(pairs_path, output_path, method="ransac")` → dict with before/after outlier counts
- **CLI:** `--input rt_pairs.tsv --output model.tsv --method ransac`
- **Tests:** TSV with 20 iRT/measured RT pairs (2 outliers), verify outlier removal and linear fit

### 10.6 `mrm_decoy_generator`
- **Wraps:** `oms.MRMDecoy()`
- **Function:** `generate_mrm_decoys(input_path, output_path, method="shuffle")` → int
- **Tests:** Small transition list, verify decoy count matches, sequences differ

**Steps:** Same as Task 1

---

## Task 11: Spectra Merging + Mass Recalibration (2 tools)

### 11.1 `spectra_merger` — `proteomics/spectrum_analysis/`
- **Wraps:** `oms.SpectraMerger()`
- **Key code:** `merger = oms.SpectraMerger(); params = merger.getDefaults(); params.setValue("mz_binning_width", 0.01); merger.setParameters(params); merger.mergeSpectraBlockWise(exp)` (for block mode) or `merger.mergeSpectraPrecursors(exp)` (for precursor mode)
- **CLI:** `--input run.mzML --output merged.mzML --mode block --block-size 3`
- **Function:** `merge_spectra(input_path, output_path, mode="block", block_size=3)` → int
- **Tests:** 9 synthetic MS1 spectra, block merge with size 3, verify 3 output spectra

### 11.2 `internal_calibrator` — `proteomics/quality_control/`
- **Wraps:** `oms.InternalCalibration()`, `oms.CalibrationData()`, `oms.MZTrafoModel()`
- **CLI:** `--input run.mzML --ids peptides.idXML --output calibrated.mzML --model linear`
- **Function:** `calibrate_mz(input_path, ids_path, output_path, model="linear")` → dict
- **Tests:** Synthetic mzML with +5 ppm systematic error, idXML with reference masses, verify error reduced after calibration

**Steps:** Same as Task 1

---

## Task 12: Search Engine + Decoy Generator (2 tools)

**Topic:** `proteomics/identification/`

### 12.1 `simple_search_engine`
- **Wraps:** `oms.SimpleSearchEngineAlgorithm()`
- **Key code:**
```python
sse = oms.SimpleSearchEngineAlgorithm()
params = sse.getDefaults()
params.setValue("precursor:mass_tolerance", precursor_tol)
params.setValue("fragment:mass_tolerance", fragment_tol)
sse.setParameters(params)
protein_ids = []
peptide_ids = []
sse.search(input_path, database_path, protein_ids, peptide_ids)
```
- **CLI:** `--input run.mzML --database proteins.fasta --output results.idXML --precursor-tol 10.0 --fragment-tol 20.0`
- **Function:** `search(input_path, database_path, output_path, precursor_tol=10.0, fragment_tol=20.0)` → int
- **Tests:** Synthetic mzML with theoretical spectrum of PEPTIDEK, small FASTA containing PEPTIDEK protein, verify identification

### 12.2 `decoy_database_generator`
- **Wraps:** `oms.DecoyGenerator()`
- **Key code:**
```python
dg = oms.DecoyGenerator()
proteins = []
oms.FASTAFile().load(input_path, proteins)
decoys = []
for entry in proteins:
    decoy = oms.FASTAEntry()
    decoy.identifier = "DECOY_" + entry.identifier
    decoy.description = entry.description
    seq = oms.AASequence.fromString(entry.sequence)
    if method == "reverse":
        rev = dg.reverseProtein(seq)
    else:
        rev = dg.shuffle(seq)
    decoy.sequence = rev.toString()
    decoys.append(decoy)
oms.FASTAFile().store(output_path, proteins + decoys)
```
- **CLI:** `--input target.fasta --output target_decoy.fasta --method reverse`
- **Function:** `generate_decoys(input_path, output_path, method="reverse")` → int
- **Tests:** 3-protein FASTA, verify output has 6 entries, decoy prefix applied, sequences reversed

**Steps:** Same as Task 1

---

## Task 13: Accurate Mass Search — Metabolomics (2 tools)

**Topic:** `metabolomics/identification/`

### 13.1 `accurate_mass_searcher`
- **Wraps:** `oms.AccurateMassSearchEngine()`
- **CLI:** `--input features.featureXML --database hmdb.tsv --output results.mzTab --mass-tol 5.0`
- **Function:** `search_accurate_mass(input_path, database_path, output_path, mass_tol=5.0)` → int
- **Tests:** Synthetic featureXML with feature at glucose [M+H]+ m/z, small TSV DB with glucose entry, verify match

### 13.2 `metabolite_spectral_matcher`
- **Wraps:** `oms.MetaboliteSpectralMatching()`
- **CLI:** `--input spectra.mzML --library lib.msp --output matches.mzTab --precursor-tol 0.1`
- **Function:** `match_spectra(input_path, library_path, output_path, precursor_tol=0.1)` → int
- **Tests:** Synthetic query spectrum matching a synthetic library spectrum, verify match score > 0

**Steps:** Same as Task 1

---

## Task 14: Metabolomics Preprocessing (2 tools)

**Topic:** `metabolomics/preprocessing/`

### 14.1 `mass_trace_detector`
- **Wraps:** `oms.MassTraceDetection()`
- **Key code:**
```python
mtd = oms.MassTraceDetection()
params = mtd.getDefaults()
params.setValue("mass_error_ppm", mass_error_ppm)
params.setValue("noise_threshold_int", noise_threshold)
mtd.setParameters(params)
mass_traces = []
mtd.run(exp, mass_traces, 0)
```
- **CLI:** `--input run.mzML --output traces.mzML --mass-error-ppm 10.0 --noise-threshold 1000`
- **Function:** `detect_mass_traces(input_path, output_path, mass_error_ppm=10.0, noise_threshold=1000.0)` → int
- **Note:** Stores mass trace data as chromatograms in output mzML
- **Tests:** Synthetic LC-MS mzML with a persistent m/z peak across multiple RT scans, verify trace detected

### 14.2 `elution_peak_detector`
- **Wraps:** `oms.MassTraceDetection()` + `oms.ElutionPeakDetection()`
- **Key code:**
```python
# Step 1: detect mass traces
mtd = oms.MassTraceDetection()
mass_traces = []
mtd.run(exp, mass_traces, 0)
# Step 2: detect elution peaks
epd = oms.ElutionPeakDetection()
params = epd.getDefaults()
params.setValue("width_filtering", width_filtering)
epd.setParameters(params)
split_traces = []
epd.detectPeaks(mass_traces, split_traces)
```
- **CLI:** `--input run.mzML --output peaks.mzML --width-filtering auto`
- **Function:** `detect_elution_peaks(input_path, output_path, width_filtering="auto")` → int
- **Note:** Runs MassTraceDetection internally then applies ElutionPeakDetection, avoiding serialization issues
- **Tests:** Synthetic LC-MS data with Gaussian elution profiles, verify peak boundaries detected

**Steps:** Same as Task 1

---

## Execution Order

Tasks 0 must run first (topic directory registration). Tasks 1-14 are fully independent and can run in parallel.

| Task | Tools | Depends On |
|------|-------|------------|
| 0 | AGENTS.md/CLAUDE.md update | — |
| 1 | gaussian_smoother, savitzky_golay_smoother, baseline_corrector, spectrum_normalizer | 0 |
| 2 | noise_filter_window, noise_filter_threshold, noise_filter_nlargest | 0 |
| 3 | signal_to_noise_estimator, peak_picker_hires, chromatogram_peak_picker, peak_type_estimator, linear_resampler | 0 |
| 4 | rt_alignment_pose_clustering, rt_alignment_identification, rt_alignment_kd | 0 |
| 5 | feature_linker_qt, feature_linker_kd, feature_linker_labeled | 0 |
| 6 | isobaric_channel_extractor, isobaric_isotope_corrector, isobaric_quantifier | 0 |
| 7 | multiplex_feature_finder, consensus_map_normalizer, id_mapper, peptide_protein_quantifier, id_conflict_resolver | 0 |
| 8 | fdr_estimator, posterior_error_probability, protein_inference_basic, protein_inference_bayesian, peptide_indexer, id_filter, consensus_id | 0 |
| 9 | deisotoper, flash_deconv, metabolite_feature_deconvolution | 0 |
| 10 | chromatogram_extractor, mrm_feature_finder, dia_scorer, mrm_transition_group_picker, mrm_rt_normalizer, mrm_decoy_generator | 0 |
| 11 | spectra_merger, internal_calibrator | 0 |
| 12 | simple_search_engine, decoy_database_generator | 0 |
| 13 | accurate_mass_searcher, metabolite_spectral_matcher | 0 |
| 14 | mass_trace_detector, elution_peak_detector | 0 |
