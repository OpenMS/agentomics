# pyopenms Gap Coverage Tools — Design Spec

## Purpose

Add 50 new standalone CLI tools that cover major pyopenms capabilities not yet represented in the agentomics tool catalog. These fill foundational gaps in signal processing, alignment, quantification, inference, DIA analysis, deconvolution, and metabolomics workflows.

All tools follow the established repo pattern: self-contained directory, click CLI, importable functions, pyopenms import guard, pytest tests with `pytest.importorskip("pyopenms")`.

### Relationship to OpenMS TOPP tools

Many of these tools wrap pyopenms classes that also have corresponding C++ OpenMS TOPP tool equivalents (e.g., `PeakPickerHiRes`, `MapAlignerPoseClustering`, `IDFilter`). This is intentional. The agentomics tools provide a **different interface**: importable Python functions with type hints for pipeline composition, click CLIs with Python-native parameter handling, and easier integration with data science workflows (pandas, notebooks). The AGENTS.md rule against duplicating TOPP tools applies to reimplementing the same algorithm in Python — not to providing a Python-native wrapper around the pyopenms bindings of that algorithm.

## Prerequisites

Before implementing the first tool in a new topic, update the topic lists in:
- `AGENTS.md` (proteomics topics and metabolomics topics sections)
- `CLAUDE.md` (topics section)

## New Topic Directories

### Proteomics

| Topic | Purpose |
|-------|---------|
| `signal_processing/` | Smoothing, baseline correction, normalization, peak picking, S/N estimation |
| `map_alignment/` | RT alignment across runs, feature linking/grouping |
| `quantification/` | TMT/iTRAQ, SILAC, label-free quantification, consensus normalization |
| `inference/` | FDR estimation, protein inference, peptide indexing, ID filtering |
| `dia_analysis/` | OpenSWATH pipeline, DIA scoring, MRM tools, chromatogram extraction |

### Metabolomics

| Topic | Purpose |
|-------|---------|
| `identification/` | Accurate mass search, spectral matching |
| `preprocessing/` | Mass trace detection, elution peak detection |

## Tool Inventory

### Return type convention

All file-I/O functions follow the existing repo pattern:
- Functions return an `int` (count of items processed) or a `dict`/`list` of summary results
- File output is a side effect (written to `output_path`)
- The CLI prints a summary line to stdout

### 1. Signal Processing (`proteomics/signal_processing/`)

#### 1.1 `gaussian_smoother`
- **Wraps:** `GaussFilter`
- **CLI:** `--input run.mzML --output smoothed.mzML --gaussian-width 0.2`
- **Functions:** `smooth_experiment(input_path, output_path, gaussian_width)` → int (spectra count)
- **Tests:** Create synthetic profile mzML with noisy peaks, verify smoothing reduces noise while preserving peak location

#### 1.2 `savitzky_golay_smoother`
- **Wraps:** `SavitzkyGolayFilter`
- **CLI:** `--input run.mzML --output smoothed.mzML --frame-length 11 --polynomial-order 3`
- **Functions:** `smooth_experiment(input_path, output_path, frame_length, polynomial_order)` → int (spectra count)
- **Tests:** Synthetic profile data, verify peak shape preservation

#### 1.3 `baseline_corrector`
- **Wraps:** `MorphologicalFilter`
- **CLI:** `--input run.mzML --output corrected.mzML --struct-element-length 3.0`
- **Functions:** `correct_baseline(input_path, output_path, struct_element_length)` → int (spectra count)
- **Tests:** Synthetic data with elevated baseline, verify baseline removal

#### 1.4 `spectrum_normalizer`
- **Wraps:** `Normalizer`
- **CLI:** `--input run.mzML --output normalized.mzML --method to_one|to_TIC`
- **Functions:** `normalize_experiment(input_path, output_path, method)` → int (spectra count)
- **Tests:** Verify max intensity = 1.0 (to_one) or sum = 1.0 (to_TIC)

#### 1.5 `noise_filter_window`
- **Wraps:** `WindowMower`
- **CLI:** `--input run.mzML --output filtered.mzML --window-size 50.0 --peak-count 3`
- **Functions:** `filter_window(input_path, output_path, window_size, peak_count)` → int (spectra count)
- **Tests:** Verify peak count reduction per window, strong peaks retained

#### 1.6 `noise_filter_threshold`
- **Wraps:** `ThresholdMower`
- **CLI:** `--input run.mzML --output filtered.mzML --threshold 100.0`
- **Functions:** `filter_threshold(input_path, output_path, threshold)` → int (spectra count)
- **Tests:** Verify all remaining peaks above threshold

#### 1.7 `noise_filter_nlargest`
- **Wraps:** `NLargest`
- **CLI:** `--input run.mzML --output filtered.mzML --n 100`
- **Functions:** `filter_nlargest(input_path, output_path, n)` → int (spectra count)
- **Tests:** Verify each spectrum has at most N peaks

#### 1.8 `signal_to_noise_estimator`
- **Wraps:** `SignalToNoiseEstimatorMedian`
- **CLI:** `--input run.mzML --output sn_report.tsv`
- **Functions:** `estimate_sn(input_path, output_path)` → int (spectra analyzed)
- **Tests:** Synthetic spectrum with known S/N, verify estimation accuracy

#### 1.9 `peak_picker_hires`
- **Wraps:** `PeakPickerHiRes`
- **CLI:** `--input profile.mzML --output centroid.mzML --signal-to-noise 1.0`
- **Functions:** `pick_peaks(input_path, output_path, signal_to_noise)` → int (spectra picked)
- **Tests:** Synthetic profile data, verify centroided output has fewer data points, correct peak positions

#### 1.10 `chromatogram_peak_picker`
- **Wraps:** `PeakPickerChromatogram`
- **CLI:** `--input run.mzML --output picked.mzML`
- **Functions:** `pick_chromatogram_peaks(input_path, output_path)` → int (chromatograms picked)
- **Tests:** Synthetic chromatogram data, verify peak detection

#### 1.11 `peak_type_estimator`
- **Wraps:** `PeakTypeEstimator`
- **CLI:** `--input run.mzML --output report.tsv`
- **Functions:** `estimate_peak_type(input_path, output_path)` → dict with type (profile/centroid) and counts
- **Tests:** Synthetic profile and centroid data, verify correct classification

#### 1.12 `linear_resampler`
- **Wraps:** `LinearResampler`
- **CLI:** `--input run.mzML --output resampled.mzML --spacing 0.01`
- **Functions:** `resample_experiment(input_path, output_path, spacing)` → int (spectra resampled)
- **Tests:** Verify uniform m/z spacing in output

### 2. Map Alignment (`proteomics/map_alignment/`)

#### 2.1 `rt_alignment_pose_clustering`
- **Wraps:** `MapAlignmentAlgorithmPoseClustering`
- **CLI:** `--reference ref.featureXML --input sample.featureXML --output aligned.featureXML --trafo trafo.transformationXML`
- **Functions:** `align_pose_clustering(reference_path, input_path, output_path)` → dict with stats (n_features, max_rt_shift)
- **Tests:** Synthetic featureXML with RT shift, verify alignment reduces RT deviation

#### 2.2 `rt_alignment_identification`
- **Wraps:** `MapAlignmentAlgorithmIdentification`
- **CLI:** `--reference ref.idXML --input sample.idXML --output aligned.idXML --model linear|b_spline|lowess`
- **Functions:** `align_identification(reference_path, input_path, output_path, model)` → dict with stats
- **Tests:** Synthetic idXML with shifted RTs, verify correction

#### 2.3 `rt_alignment_kd`
- **Wraps:** `MapAlignmentAlgorithmKD`
- **CLI:** `--inputs run1.featureXML run2.featureXML --output-dir aligned/ --model lowess`
- **Functions:** `align_kd(input_paths, output_dir, model)` → dict with stats per file
- **Tests:** Multiple synthetic featureXMLs with offsets, verify convergence

#### 2.4 `feature_linker_qt`
- **Wraps:** `FeatureGroupingAlgorithmQT`
- **CLI:** `--inputs run1.featureXML run2.featureXML --output consensus.consensusXML --rt-tol 30.0 --mz-tol 10.0`
- **Functions:** `link_features_qt(input_paths, output_path, rt_tol, mz_tol)` → int (consensus features)
- **Tests:** Synthetic featureXMLs with shared features, verify correct grouping

#### 2.5 `feature_linker_kd`
- **Wraps:** `FeatureGroupingAlgorithmKD`
- **CLI:** `--inputs run1.featureXML run2.featureXML --output consensus.consensusXML --rt-tol 30.0 --mz-tol 10.0`
- **Functions:** `link_features_kd(input_paths, output_path, rt_tol, mz_tol)` → int (consensus features)
- **Tests:** Same as QT linker tests

#### 2.6 `feature_linker_labeled`
- **Wraps:** `FeatureGroupingAlgorithmLabeled`
- **CLI:** `--input run.featureXML --output consensus.consensusXML --rt-estimate 1000.0 --mz-estimate 5.0`
- **Functions:** `link_labeled_features(input_path, output_path, rt_estimate, mz_estimate)` → int (consensus features)
- **Tests:** Synthetic labeled feature pairs (light/heavy), verify pairing

### 3. Quantification (`proteomics/quantification/`)

**Isobaric workflow order:** `isobaric_channel_extractor` → `isobaric_isotope_corrector` → `isobaric_quantifier`. Each tool's output feeds the next.

#### 3.1 `isobaric_channel_extractor`
- **Wraps:** `IsobaricChannelExtractor`
- **CLI:** `--input run.mzML --method tmt6plex|tmt10plex|tmt16plex|itraq4plex|itraq8plex --output quant.consensusXML`
- **Functions:** `extract_channels(input_path, method, output_path)` → int (spectra with reporter ions)
- **Tests:** Synthetic MS2/MS3 with reporter ions at known m/z, verify extraction

#### 3.2 `isobaric_isotope_corrector`
- **Wraps:** `IsobaricIsotopeCorrector`
- **CLI:** `--input quant.consensusXML --method tmt6plex --output corrected.consensusXML`
- **Functions:** `correct_isotope_impurities(input_path, method, output_path)` → int (features corrected)
- **Tests:** Known cross-channel contamination, verify correction

#### 3.3 `isobaric_quantifier`
- **Wraps:** `IsobaricQuantifier`
- **CLI:** `--input quant.consensusXML --output quantified.consensusXML`
- **Functions:** `quantify_isobaric(input_path, output_path)` → int (quantified features)
- **Tests:** Synthetic consensus data, verify quantification output

#### 3.4 `multiplex_feature_finder`
- **Wraps:** `FeatureFinderMultiplexAlgorithm`
- **CLI:** `--input run.mzML --output features.featureXML --labels light:0 heavy:8.0142`
- **Functions:** `find_multiplex_features(input_path, output_path, labels)` → int (feature pairs/triplets found)
- **Tests:** Synthetic SILAC-like data with known mass shifts, verify pair detection

#### 3.5 `consensus_map_normalizer`
- **Wraps:** `ConsensusMapNormalizerAlgorithmMedian`, `ConsensusMapNormalizerAlgorithmQuantile`
- **CLI:** `--input consensus.consensusXML --output normalized.consensusXML --method median|quantile`
- **Functions:** `normalize_consensus(input_path, output_path, method)` → int (features normalized)
- **Tests:** Synthetic consensus map with biased intensities, verify normalization

#### 3.6 `id_mapper`
- **Wraps:** `IDMapper`
- **CLI:** `--features features.featureXML --ids peptides.idXML --output annotated.featureXML --rt-tol 5.0 --mz-tol 10.0`
- **Functions:** `map_ids_to_features(features_path, ids_path, output_path, rt_tol, mz_tol)` → int (IDs mapped)
- **Tests:** Synthetic features and IDs with matching RT/mz, verify annotation

#### 3.7 `peptide_protein_quantifier`
- **Wraps:** `PeptideAndProteinQuant`
- **CLI:** `--input annotated.consensusXML --output protein_quant.csv`
- **Functions:** `quantify_proteins(input_path, output_path)` → int (proteins quantified)
- **Tests:** Synthetic consensus map with peptide annotations, verify protein rollup

#### 3.8 `id_conflict_resolver`
- **Wraps:** `IDConflictResolverAlgorithm`
- **CLI:** `--input features.featureXML --output resolved.featureXML`
- **Functions:** `resolve_id_conflicts(input_path, output_path)` → int (conflicts resolved)
- **Tests:** Features with multiple conflicting IDs, verify single best ID retained

### 4. Inference (`proteomics/inference/`)

#### 4.1 `fdr_estimator`
- **Wraps:** `FalseDiscoveryRate`
- **CLI:** `--input peptides.idXML --output fdr.idXML --protein/--no-protein`
- **Functions:** `estimate_fdr(input_path, output_path, protein_level)` → dict with counts at FDR thresholds
- **Tests:** Synthetic idXML with target/decoy hits, verify q-value calculation

#### 4.2 `posterior_error_probability`
- **Wraps:** `PosteriorErrorProbabilityModel`
- **CLI:** `--input peptides.idXML --output pep.idXML`
- **Functions:** `estimate_pep(input_path, output_path)` → int (PSMs scored)
- **Tests:** Synthetic score distributions, verify PEP assignment

#### 4.3 `protein_inference_basic`
- **Wraps:** `BasicProteinInferenceAlgorithm`
- **CLI:** `--input peptides.idXML --output proteins.idXML`
- **Functions:** `infer_proteins(input_path, output_path)` → int (proteins inferred)
- **Tests:** Peptides mapping to known proteins, verify correct inference

#### 4.4 `protein_inference_bayesian`
- **Wraps:** `BayesianProteinInferenceAlgorithm`
- **CLI:** `--input peptides.idXML --output proteins.idXML`
- **Functions:** `infer_proteins_bayesian(input_path, output_path)` → int (proteins inferred)
- **Tests:** Shared peptide scenario, verify Bayesian handling

#### 4.5 `peptide_indexer`
- **Wraps:** `PeptideIndexing`
- **CLI:** `--ids peptides.idXML --fasta database.fasta --output indexed.idXML`
- **Functions:** `index_peptides(ids_path, fasta_path, output_path)` → int (peptides indexed)
- **Tests:** Known peptides from known proteins, verify correct mapping

#### 4.6 `id_filter`
- **Wraps:** `IDFilter`
- **CLI:** `--input peptides.idXML --output filtered.idXML --score-threshold 0.05 --remove-decoys`
- **Functions:** `filter_ids(input_path, output_path, score_threshold, remove_decoys)` → dict with before/after counts
- **Tests:** Mixed target/decoy with varied scores, verify correct filtering

#### 4.7 `consensus_id`
- **Wraps:** `ConsensusIDAlgorithmAverage`, `ConsensusIDAlgorithmBest`, `ConsensusIDAlgorithmRanks`
- **CLI:** `--inputs search1.idXML search2.idXML --output consensus.idXML --algorithm average|best|ranks`
- **Functions:** `consensus_id(input_paths, output_path, algorithm)` → int (consensus PSMs)
- **Tests:** Two search results with overlapping IDs, verify consensus

### 5. Charge/Isotope Deconvolution (`proteomics/spectrum_analysis/`)

These fit the existing `spectrum_analysis/` topic.

#### 5.1 `deisotoper`
- **Wraps:** `Deisotoper`
- **CLI:** `--input run.mzML --output deisotoped.mzML --fragment-tolerance 0.1 --min-charge 1 --max-charge 5`
- **Functions:** `deisotope(input_path, output_path, fragment_tolerance, min_charge, max_charge)` → int (spectra processed)
- **Tests:** Synthetic spectrum with isotope envelopes, verify monoisotopic reduction

#### 5.2 `flash_deconv`
- **Wraps:** `FLASHDeconvAlgorithm`
- **CLI:** `--input intact.mzML --output masses.tsv --min-mass 5000 --max-mass 100000`
- **Functions:** `deconvolve_intact(input_path, output_path, min_mass, max_mass)` → int (masses found)
- **Note:** Higher complexity implementation — FLASHDeconv has many parameters and generating realistic intact protein test data is non-trivial
- **Tests:** Synthetic multiply-charged protein spectrum, verify mass recovery

### 6. Metabolite Feature Deconvolution (`metabolomics/feature_processing/`)

Fits the existing `feature_processing/` topic since it groups adducts in metabolite feature maps.

#### 6.1 `metabolite_feature_deconvolution`
- **Wraps:** `MetaboliteFeatureDeconvolution`
- **CLI:** `--input features.featureXML --output grouped.consensusXML`
- **Functions:** `deconvolve_adducts(input_path, output_path)` → int (adduct groups found)
- **Tests:** Synthetic features at [M+H]+, [M+Na]+, [M+K]+, verify grouping

### 7. OpenSWATH / DIA (`proteomics/dia_analysis/`)

Separate from `targeted_proteomics/` because these are OpenSWATH pipeline-specific algorithms (scoring, decoy generation, transition group picking) rather than general targeted/SRM utilities.

#### 7.1 `chromatogram_extractor`
- **Wraps:** `ChromatogramExtractor`, `ChromatogramExtractorAlgorithm`
- **CLI:** `--input dia.mzML --transitions transitions.tsv --output chromatograms.mzML --mz-tol 10.0`
- **Functions:** `extract_chromatograms(input_path, transitions_path, output_path, mz_tol)` → int (chromatograms extracted)
- **Tests:** Synthetic DIA data with known transitions, verify XIC extraction

#### 7.2 `mrm_feature_finder`
- **Wraps:** `MRMFeatureFinderScoring`
- **CLI:** `--input chromatograms.mzML --transitions transitions.tsv --output features.featureXML`
- **Functions:** `find_mrm_features(input_path, transitions_path, output_path)` → int (feature groups found)
- **Tests:** Synthetic chromatogram with clear peaks, verify feature detection

#### 7.3 `dia_scorer`
- **Wraps:** `DIAScoring`
- **CLI:** `--input dia.mzML --transitions transitions.tsv --output scores.tsv`
- **Functions:** `score_dia(input_path, transitions_path, output_path)` → int (transitions scored)
- **Tests:** Synthetic DIA spectrum, verify score computation

#### 7.4 `mrm_transition_group_picker`
- **Wraps:** `MRMTransitionGroupPicker`
- **CLI:** `--input chromatograms.mzML --output picked.mzML`
- **Functions:** `pick_transition_groups(input_path, output_path)` → int (groups picked)
- **Tests:** Synthetic transition group chromatograms, verify peak picking

#### 7.5 `mrm_rt_normalizer`
- **Wraps:** `MRMRTNormalizer`
- **CLI:** `--input features.featureXML --output normalized.featureXML`
- **Functions:** `normalize_mrm_rt(input_path, output_path)` → int (features normalized)
- **Tests:** Features with known iRT peptides, verify normalization

#### 7.6 `mrm_decoy_generator`
- **Wraps:** `MRMDecoy`
- **CLI:** `--input transitions.tsv --output decoy_transitions.tsv --method shuffle|reverse|pseudo-reverse`
- **Functions:** `generate_mrm_decoys(input_path, output_path, method)` → int (decoy transitions generated)
- **Tests:** Known transition list, verify decoy generation properties (different sequences, same precursor properties)

### 8. Spectra Merging (`proteomics/spectrum_analysis/`)

#### 8.1 `spectra_merger`
- **Wraps:** `SpectraMerger`
- **CLI:** `--input run.mzML --output merged.mzML --mode block|precursor --block-size 3`
- **Functions:** `merge_spectra(input_path, output_path, mode, block_size)` → int (spectra after merging)
- **Tests:** Synthetic consecutive spectra, verify merging reduces spectrum count, preserves peaks

### 9. Mass Recalibration (`proteomics/quality_control/`)

Fits the existing `quality_control/` topic.

#### 9.1 `internal_calibrator`
- **Wraps:** `InternalCalibration`, `CalibrationData`, `MZTrafoModel`
- **CLI:** `--input run.mzML --ids peptides.idXML --output calibrated.mzML --model linear|quadratic`
- **Functions:** `calibrate_mz(input_path, ids_path, output_path, model)` → dict with before/after ppm stats
- **Tests:** Synthetic mzML with systematic m/z error, synthetic idXML with reference masses, verify error reduction

### 10. Search Engine (`proteomics/identification/`)

Fits the existing `identification/` topic.

#### 10.1 `simple_search_engine`
- **Wraps:** `SimpleSearchEngineAlgorithm`
- **CLI:** `--input run.mzML --database proteins.fasta --output results.idXML --precursor-tol 10.0 --fragment-tol 20.0`
- **Functions:** `search(input_path, database_path, output_path, precursor_tol, fragment_tol)` → int (PSMs found)
- **Tests:** Synthetic mzML with known peptide spectra, small FASTA, verify correct identification

#### 10.2 `decoy_database_generator`
- **Wraps:** `DecoyGenerator`
- **CLI:** `--input target.fasta --output target_decoy.fasta --method reverse|shuffle`
- **Functions:** `generate_decoys(input_path, output_path, method)` → int (decoy proteins generated)
- **Tests:** Small FASTA, verify decoy count equals target, prefix applied, sequences differ

### 11. Accurate Mass Search — Metabolomics (`metabolomics/identification/`)

#### 11.1 `accurate_mass_searcher`
- **Wraps:** `AccurateMassSearchEngine`
- **CLI:** `--input features.featureXML --database hmdb.tsv --output results.mzTab --mass-tol 5.0`
- **Functions:** `search_accurate_mass(input_path, database_path, output_path, mass_tol)` → int (matches found)
- **Tests:** Synthetic features at known metabolite masses, small DB, verify matches

#### 11.2 `metabolite_spectral_matcher`
- **Wraps:** `MetaboliteSpectralMatching`
- **CLI:** `--input spectra.mzML --library lib.msp --output matches.mzTab --precursor-tol 0.1`
- **Functions:** `match_spectra(input_path, library_path, output_path, precursor_tol)` → int (spectra matched)
- **Tests:** Synthetic query spectrum matching library entry, verify match score

### 12. Metabolomics Preprocessing (`metabolomics/preprocessing/`)

**Relationship to existing tools:** The existing `metabolite_feature_detection` tool in `metabolomics/feature_processing/` runs `MassTraceDetection` + `ElutionPeakDetection` as a monolithic pipeline. These tools decompose that into individually usable steps for finer-grained control.

#### 12.1 `mass_trace_detector`
- **Wraps:** `MassTraceDetection`
- **CLI:** `--input run.mzML --output traces.mzML --mass-error-ppm 10.0 --noise-threshold 1000`
- **Functions:** `detect_mass_traces(input_path, output_path, mass_error_ppm, noise_threshold)` → int (traces found)
- **Note:** Outputs an mzML with mass traces stored as chromatograms, not TSV, since `MassTrace` objects are not trivially serializable to flat files
- **Tests:** Synthetic LC-MS data with known compounds, verify trace detection at expected m/z and RT

#### 12.2 `elution_peak_detector`
- **Wraps:** `ElutionPeakDetection`
- **CLI:** `--input run.mzML --output peaks.mzML --width-filtering auto`
- **Functions:** `detect_elution_peaks(input_path, output_path, width_filtering)` → int (peaks found)
- **Note:** Takes mzML as input, runs `MassTraceDetection` internally to produce `MassTrace` objects, then applies `ElutionPeakDetection`. This avoids the serialization problem while still being independently useful.
- **Tests:** Synthetic LC-MS data with Gaussian elution profiles, verify peak boundaries

## Summary

| Area | Topic | Tool Count |
|------|-------|-----------|
| Signal Processing | `proteomics/signal_processing/` | 12 |
| Map Alignment + Feature Linking | `proteomics/map_alignment/` | 6 |
| Quantification | `proteomics/quantification/` | 8 |
| Inference / FDR | `proteomics/inference/` | 7 |
| Deconvolution | `proteomics/spectrum_analysis/` | 2 |
| Metabolite Deconvolution | `metabolomics/feature_processing/` | 1 |
| OpenSWATH / DIA | `proteomics/dia_analysis/` | 6 |
| Spectra Merging | `proteomics/spectrum_analysis/` | 1 |
| Mass Recalibration | `proteomics/quality_control/` | 1 |
| Search Engine | `proteomics/identification/` | 2 |
| Accurate Mass Search | `metabolomics/identification/` | 2 |
| Metabolomics Preprocessing | `metabolomics/preprocessing/` | 2 |
| **Total** | | **50** |

## Notes

- `FeatureFinderMetabo` is not available in the installed pyopenms — the metabolomics preprocessing tools will use `MassTraceDetection` + `ElutionPeakDetection` directly instead
- All 79 other pyopenms classes verified as available
- Tools that fit existing topics (deconvolution → `spectrum_analysis/`, recalibration → `quality_control/`, search → `identification/`, metabolite deconvolution → `feature_processing/`) go there rather than creating new topics
- The `noise_filter` was split into 3 separate tools (`noise_filter_window`, `noise_filter_threshold`, `noise_filter_nlargest`) because `WindowMower`, `ThresholdMower`, and `NLargest` have different parameter semantics
