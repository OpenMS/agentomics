# Agentomics Tool Catalog — 145 CLI Tools for Proteomics & Metabolomics

**Date:** 2026-03-24
**Status:** Draft — awaiting user approval
**Scope:** Only tools providing NEW functionality not available as a simple OpenMS TOPP command
**Tools 1-100:** Core utilities derived from GitHub popularity, community requests, pyopenms docs
**Tools 101-145:** Niche paper-level tools from published research workflows

## Scoring Criteria

Each tool scored 1-10 on:
- **Popularity** (P): GitHub stars/forks, community requests, forum frequency
- **Uniqueness** (U): Not available as a TOPP tool or simple pyopenms one-liner
- **Utility** (V): How often researchers need this in daily workflows
- **Final Score** = (P + U + V) / 3, rounded

---

## Already Implemented (7 tools — excluded from catalog)

| # | Tool | Domain |
|---|------|--------|
| - | peptide_mass_calculator | proteomics |
| - | protein_digest | proteomics |
| - | feature_detection_proteomics | proteomics |
| - | spectrum_file_info | proteomics |
| - | isotope_pattern_matcher | metabolomics |
| - | mass_accuracy_calculator | metabolomics |
| - | metabolite_feature_detection | metabolomics |

---

## Excluded — Direct TOPP Tool Duplicates

These exist as standalone TOPP commands and are excluded per option A:

FileConverter, PeakPickerHiRes, NoiseFilterGaussian, NoiseFilterSGolay, Normalizer, SpectraMerger, DecoyDatabase, SimpleSearchEngine, FalseDiscoveryRate, IDFilter, IDMapper, ProteinInference, PeptideIndexing, MapAlignerPoseClustering, MapAlignerIdentification, FeatureLinkerUnlabeledKD, FeatureLinkerUnlabeledQT, ConsensusMapNormalizer, IsobaricAnalyzer, InternalCalibration, GNPSExport, AccurateMassSearch, MetaboliteAdductDecharger, FeatureFinderMetabo, FeatureFinderCentroided, IDFileConverter, TextExporter, MRMMapper, OpenSwathChromatogramExtractor, OpenSwathRTNormalizer, QCCalculator, FileFilter, BaselineFilter

---

## Category 1: Spectrum Analysis & Annotation (15 tools)

### 1. theoretical_spectrum_generator
- **Description:** Generate theoretical b/y/a/c/x/z fragment ion spectra for a peptide sequence with annotated ion labels. Outputs human-readable TSV with ion type, number, charge, m/z, and annotation. Unlike any TOPP tool — TOPP has no standalone fragment spectrum output.
- **pyopenms classes:** `TheoreticalSpectrumGenerator`, `AASequence`, `MSSpectrum`, `Param`
- **CLI:** `--sequence PEPTM(Oxidation)IDEK --charge 2 --ion-types b,y,a --add-losses --add-isotopes --output fragments.tsv`
- **Inputs:** Peptide sequence (ProForma/bracket notation), charge state, ion type selection
- **Outputs:** TSV (ion_type, ion_number, charge, mz, annotation, intensity)
- **Domain:** proteomics
- **Why needed:** Universal need for PSM validation, teaching, SRM/PRM method development, spectral library QC. spectrum_utils (164 stars) addresses this but requires heavy dependencies.
- **Score:** P=8, U=8, V=9 → **8**

### 2. spectrum_similarity_scorer
- **Description:** Compute pairwise spectral similarity (cosine, modified cosine, spectral contrast angle) between MS2 spectra from mzML or MGF files. No TOPP tool computes and outputs similarity scores.
- **pyopenms classes:** `SpectrumAlignment`, `MSSpectrum`, `MSExperiment`, `MzMLFile`
- **CLI:** `--query query.mgf --library reference.mgf --tolerance 0.02 --metric cosine --output scores.tsv`
- **Inputs:** Two spectrum files (mzML or MGF), tolerance, metric choice
- **Outputs:** TSV (query_id, library_id, score, matched_peaks, total_query_peaks, total_lib_peaks)
- **Domain:** proteomics / metabolomics
- **Why needed:** matchms (253 stars, 77 forks) is the go-to but requires a heavy install. A pyopenms-native CLI fills this gap. Fundamental for molecular networking and library matching.
- **Score:** P=9, U=9, V=9 → **9**

### 3. spectrum_annotator
- **Description:** Given an observed MS2 spectrum and a peptide sequence, annotate peaks with theoretical fragment ion matches. Output annotation data for mirror plots or publication figures.
- **pyopenms classes:** `TheoreticalSpectrumGenerator`, `SpectrumAlignment`, `AASequence`, `MSSpectrum`
- **CLI:** `--spectrum observed.mzML --scan-index 1234 --sequence PEPTIDEK --charge 2 --tolerance 0.02 --output annotation.tsv`
- **Inputs:** mzML + scan index (or MGF + scan), peptide sequence, charge, tolerance
- **Outputs:** TSV (mz, intensity, annotation, matched, error_da, error_ppm)
- **Domain:** proteomics
- **Why needed:** spectrum_utils (164 stars) does this but no simple pyopenms CLI exists. Researchers need annotated spectra for publications and QC.
- **Score:** P=8, U=8, V=8 → **8**

### 4. spectrum_scoring_hyperscore
- **Description:** Score an experimental MS2 spectrum against a theoretical spectrum using HyperScore. Returns score, matched ions, and sequence coverage.
- **pyopenms classes:** `HyperScore`, `TheoreticalSpectrumGenerator`, `AASequence`, `MSSpectrum`
- **CLI:** `--input spectrum.mzML --scan-index 100 --sequence PEPTIDEK --charge 2 --tolerance 0.02 --output score.json`
- **Inputs:** mzML + scan, peptide sequence, charge
- **Outputs:** JSON (hyperscore, matched_b, matched_y, total_matched, sequence_coverage)
- **Domain:** proteomics
- **Why needed:** Quick PSM validation without running a full search engine. No TOPP tool outputs HyperScore for a single spectrum-peptide pair.
- **Score:** P=6, U=8, V=7 → **7**

### 5. neutral_loss_scanner
- **Description:** Scan MS2 spectra for characteristic neutral losses (e.g., -98 Da phospho, -162 Da hexose, -176 Da glucuronide). Reports precursor-fragment pairs matching user-specified losses.
- **pyopenms classes:** `MSExperiment`, `MzMLFile`, `MSSpectrum`
- **CLI:** `--input file.mzML --losses 97.977,162.053,176.032 --tolerance 0.02 --output matches.tsv`
- **Inputs:** mzML file, list of neutral loss masses, tolerance
- **Outputs:** TSV (scan_id, precursor_mz, fragment_mz, neutral_loss, loss_name, intensity)
- **Domain:** metabolomics / proteomics (PTMs)
- **Why needed:** Standard metabolomics annotation strategy for drug metabolites and PTM screening. No TOPP tool or simple CLI exists.
- **Score:** P=7, U=9, V=8 → **8**

### 6. mass_defect_filter
- **Description:** Compute mass defect (fractional mass) for features, filter by mass defect range, and generate Kendrick mass defect data for homologous series detection.
- **pyopenms classes:** `EmpiricalFormula`, `MSExperiment`, `MSSpectrum`
- **CLI:** `--input features.tsv --mdf-min 0.1 --mdf-max 0.3 --kendrick-base CH2 --output filtered.tsv`
- **Inputs:** Feature list (TSV with m/z column) or mzML, mass defect range, Kendrick base
- **Outputs:** TSV (mz, mass_defect, kendrick_mass, kendrick_mass_defect, passes_filter)
- **Domain:** metabolomics
- **Why needed:** Key technique in drug metabolism and environmental metabolomics. Mass-Suite paper demonstrates need. No TOPP tool.
- **Score:** P=6, U=9, V=7 → **7**

### 7. mass_decomposition_tool
- **Description:** Given an accurate mass, find all possible molecular formula or amino acid compositions within tolerance using the MassDecompositionAlgorithm.
- **pyopenms classes:** `MassDecompositionAlgorithm`, `MassDecomposition`, `Param`
- **CLI:** `--mass 1234.567 --tolerance 0.01 --residue-set full --output decompositions.tsv`
- **Inputs:** Target mass, tolerance, element/residue constraints
- **Outputs:** TSV of candidate compositions with exact mass and error
- **Domain:** metabolomics / proteomics
- **Why needed:** Unique algorithmic capability. pyopenms docs show it as a CLI example. No TOPP tool wraps this.
- **Score:** P=7, U=9, V=7 → **8**

### 8. molecular_formula_finder
- **Description:** Enumerate valid molecular formulas for an accurate mass within ppm tolerance, constrained by element ranges. Apply Seven Golden Rules filtering (DBE, element ratios).
- **pyopenms classes:** `EmpiricalFormula`, `Element`, `IsotopeDistribution`
- **CLI:** `--mass 180.0634 --ppm 5 --elements C:0-12,H:0-30,N:0-5,O:0-10 --rules --output formulas.tsv`
- **Inputs:** Target mass, ppm tolerance, element ranges, optional golden rules filter
- **Outputs:** TSV (formula, exact_mass, error_ppm, dbe, h_c_ratio, isotope_pattern)
- **Domain:** metabolomics
- **Why needed:** Fundamental metabolomics annotation step. MetaboAnalyst and SIRIUS address this but no simple offline CLI. No TOPP tool.
- **Score:** P=8, U=8, V=8 → **8**

### 9. xic_extractor
- **Description:** Extract ion chromatograms (XIC/EIC) for target m/z values with ppm tolerance from mzML files. Also compute TIC and BPC. Outputs RT vs intensity traces.
- **pyopenms classes:** `MSExperiment`, `MzMLFile`, `MSSpectrum`
- **CLI:** `--input run.mzML --mz 524.265 --ppm 10 --rt-start 10 --rt-end 60 --output xic.tsv` or `--mz-list targets.tsv --output xics.tsv`
- **Inputs:** mzML file, target m/z (or list), ppm tolerance, optional RT range
- **Outputs:** TSV (rt, intensity) per target; optional TIC/BPC
- **Domain:** proteomics / metabolomics
- **Why needed:** One of the most common ad-hoc MS tasks. OpenMS issue #4977 documents the gap. OpenSwathChromatogramExtractor is DIA-specific and requires TraML — this is simpler.
- **Score:** P=9, U=8, V=10 → **9**

### 10. tic_bpc_calculator
- **Description:** Compute Total Ion Current and Base Peak Chromatogram from mzML, export as TSV. Separate from spectrum_file_info which reports summary stats, not full chromatogram traces.
- **pyopenms classes:** `MSExperiment`, `MzMLFile`, `MSSpectrum`
- **CLI:** `--input run.mzML --ms-level 1 --output chromatograms.tsv`
- **Inputs:** mzML file, MS level selection
- **Outputs:** TSV (rt, tic_intensity, bpc_mz, bpc_intensity)
- **Domain:** proteomics / metabolomics
- **Why needed:** Foundation for QC visualization. No TOPP tool outputs full TIC/BPC traces as flat files.
- **Score:** P=7, U=7, V=8 → **7**

### 11. dia_window_analyzer
- **Description:** Read DIA mzML and report isolation window scheme: window widths, centers, overlaps, cycle time, number of windows per cycle.
- **pyopenms classes:** `MSExperiment`, `MSSpectrum`, `Precursor`
- **CLI:** `--input dia_run.mzML --output windows.tsv`
- **Inputs:** DIA mzML file
- **Outputs:** TSV (window_center, lower_offset, upper_offset, width, overlap_with_next), summary stats
- **Domain:** proteomics
- **Why needed:** DIA analysis tools require window schemes but they're not always documented. Common pain point when setting up DIA-NN or OpenSWATH. No TOPP tool.
- **Score:** P=7, U=9, V=7 → **8**

### 12. precursor_charge_distribution
- **Description:** Analyze precursor charge state distribution across all MS2 spectra in an mzML file. Report counts and percentages per charge state.
- **pyopenms classes:** `MSExperiment`, `MzMLFile`, `MSSpectrum`, `Precursor`
- **CLI:** `--input run.mzML --output charge_dist.tsv`
- **Inputs:** mzML file
- **Outputs:** TSV (charge_state, count, percentage)
- **Domain:** proteomics
- **Why needed:** Key QC metric for evaluating ionization and instrument performance. Not directly available from any TOPP tool as standalone output.
- **Score:** P=6, U=8, V=7 → **7**

### 13. mzml_spectrum_subsetter
- **Description:** Extract specific spectra from mzML by scan number list, RT list, or precursor m/z list. Unlike FileFilter (which filters by ranges), this extracts specific spectra by ID.
- **pyopenms classes:** `MSExperiment`, `MzMLFile`, `MSSpectrum`
- **CLI:** `--input run.mzML --scans 100,200,300 --output subset.mzML` or `--precursor-mz-list targets.tsv --output subset.mzML`
- **Inputs:** mzML + scan numbers or precursor m/z list
- **Outputs:** mzML with only selected spectra
- **Domain:** proteomics / metabolomics
- **Why needed:** FileFilter works on ranges; this works on specific lists. Common need for extracting spectra around identified PSMs for visualization.
- **Score:** P=6, U=7, V=7 → **7**

### 14. spectrum_entropy_calculator
- **Description:** Calculate spectral entropy and normalized entropy for MS2 spectra. Spectral entropy is emerging as a quality metric superior to simple peak count.
- **pyopenms classes:** `MSExperiment`, `MSSpectrum`, `MzMLFile`
- **CLI:** `--input run.mzML --ms-level 2 --output entropy.tsv`
- **Inputs:** mzML file
- **Outputs:** TSV (scan_index, rt, precursor_mz, entropy, normalized_entropy, n_peaks)
- **Domain:** proteomics / metabolomics
- **Why needed:** Spectral entropy is a powerful QC and scoring metric (Li et al., Nature Methods 2021). No TOPP tool computes it.
- **Score:** P=7, U=9, V=7 → **8**

### 15. spectral_library_format_converter
- **Description:** Convert between spectral library formats: MSP, sptxt, TraML, PQP/TSV (OpenSWATH). No single TOPP tool handles all conversions.
- **pyopenms classes:** `TraMLFile`, `TargetedExperiment`, `TransitionPQPFile`
- **CLI:** `--input library.msp --output library.traml --format traml`
- **Inputs:** Spectral library file (MSP, sptxt, TraML, PQP)
- **Outputs:** Converted library in target format
- **Domain:** proteomics
- **Why needed:** DIA workflows require libraries in specific formats. Converting between them is a constant pain point. No single TOPP tool does all conversions.
- **Score:** P=7, U=8, V=8 → **8**

---

## Category 2: Peptide & Protein Analysis (16 tools)

### 16. peptide_property_calculator
- **Description:** Calculate physicochemical properties: pI, hydrophobicity (multiple scales: Kyte-Doolittle, GRAVY, Hopp-Woods), charge at pH, instability index, amino acid composition.
- **pyopenms classes:** `AASequence`, `Residue`, `EmpiricalFormula`
- **CLI:** `--sequence PEPTIDEK --ph 7.0 --output properties.json` or `--input peptides.tsv --output properties.tsv`
- **Inputs:** Peptide sequence(s), pH for charge calculation
- **Outputs:** TSV/JSON (sequence, mw, pI, gravy, charge_at_ph, instability_index, aa_composition)
- **Domain:** proteomics
- **Why needed:** modlAMP (61 stars) and PepFun address this but require separate installs. No pyopenms/TOPP tool computes these properties. Essential for targeted proteomics assay design.
- **Score:** P=8, U=9, V=9 → **9**

### 17. peptide_uniqueness_checker
- **Description:** Check which peptides in a list are proteotypic (unique to a single protein) within a FASTA database.
- **pyopenms classes:** `FASTAFile`, `AASequence`, `ProteaseDigestion`
- **CLI:** `--peptides peptide_list.tsv --fasta uniprot.fasta --output uniqueness.tsv`
- **Inputs:** Peptide list (TSV), FASTA database
- **Outputs:** TSV (peptide, protein_accessions, is_unique, n_proteins)
- **Domain:** proteomics
- **Why needed:** neXtProt checker (Oxford Academic paper, 2017) is web-only. Essential for SRM/PRM assay design and HPP protein evidence. No TOPP tool.
- **Score:** P=8, U=9, V=8 → **8**

### 18. protein_coverage_calculator
- **Description:** Map identified peptides to protein sequences and calculate per-protein sequence coverage, with position mapping.
- **pyopenms classes:** `FASTAFile`, `AASequence`, `IdXMLFile` (optional)
- **CLI:** `--fasta proteins.fasta --peptides identified.tsv --accession P12345 --output coverage.tsv`
- **Inputs:** FASTA database, peptide list (TSV or idXML)
- **Outputs:** TSV (accession, length, covered_residues, coverage_pct, peptide_count, coverage_map)
- **Domain:** proteomics
- **Why needed:** PrIntMap-R exists but is R/web-only. AlphaMap (92 stars) is GUI-focused. No simple Python CLI. No TOPP tool.
- **Score:** P=8, U=9, V=8 → **8**

### 19. modification_mass_calculator
- **Description:** Query Unimod/PSI-MOD databases by name or mass shift, compute modified peptide masses, list all known modifications for a residue.
- **pyopenms classes:** `ModificationsDB`, `ResidueModification`, `AASequence`
- **CLI:** `--search-mod Phospho --list-mods` or `--sequence PEPTIDEK --modifications "Oxidation(M):4,Phospho(S):2" --charge 2`
- **Inputs:** Modification name or mass shift query, or peptide with modifications
- **Outputs:** Modification details (name, mass_shift, formula_delta, sites, classification) or modified peptide mass/mz
- **Domain:** proteomics
- **Why needed:** Researchers search Unimod manually daily. MSModDetector paper (2024) shows the need. No TOPP tool for mod database queries.
- **Score:** P=8, U=9, V=8 → **8**

### 20. modified_peptide_generator
- **Description:** Generate all modified peptide variants for given variable/fixed modifications, respecting max modifications per peptide.
- **pyopenms classes:** `ModifiedPeptideGenerator`, `AASequence`, `ModificationsDB`
- **CLI:** `--sequence PEPTMIDEK --variable-mods Oxidation,Phospho --fixed-mods Carbamidomethyl --max-mods 2 --output variants.tsv`
- **Inputs:** Peptide sequence, variable mods, fixed mods, max modifications
- **Outputs:** TSV (modified_sequence, mono_mass, mz_at_charge, n_mods, mod_positions)
- **Domain:** proteomics
- **Why needed:** Needed for building targeted method inclusion lists and search space estimation. No TOPP tool enumerates variants.
- **Score:** P=6, U=9, V=7 → **7**

### 21. peptide_modification_analyzer
- **Description:** Parse a modified peptide sequence (bracket/ProForma notation) and output residue-by-residue mass breakdown with modification details.
- **pyopenms classes:** `AASequence`, `ModificationsDB`, `ResidueDB`
- **CLI:** `--sequence ".(Dimethyl)PEPTM(Oxidation)IDE." --charge 2`
- **Inputs:** Modified peptide string, charge
- **Outputs:** JSON/TSV (residue, position, base_mass, mod_name, mod_mass, cumulative_b_ion, cumulative_y_ion)
- **Domain:** proteomics
- **Why needed:** Useful for understanding and validating search engine notation. No TOPP tool provides residue-level breakdown.
- **Score:** P=6, U=8, V=7 → **7**

### 22. missed_cleavage_analyzer
- **Description:** Analyze missed cleavage distribution in identification results. Key QC metric for digestion efficiency.
- **pyopenms classes:** `ProteaseDigestion`, `AASequence`, `PeptideIdentification` (optional)
- **CLI:** `--input results.idXML --enzyme Trypsin --output mc_report.tsv` or `--peptides peptide_list.tsv --enzyme Trypsin`
- **Inputs:** Identification results (idXML or peptide list TSV), enzyme name
- **Outputs:** TSV (missed_cleavages, count, percentage), summary stats
- **Domain:** proteomics
- **Why needed:** Standard QC metric in every proteomics facility. Calculated manually or in R. No TOPP tool. Papers cite its importance for quantification accuracy.
- **Score:** P=7, U=9, V=8 → **8**

### 23. peptide_detectability_predictor
- **Description:** Given a FASTA and digestion parameters, predict peptide detectability based on physicochemical heuristics (length, hydrophobicity, charge, MW range).
- **pyopenms classes:** `ProteaseDigestion`, `AASequence`, `FASTAFile`, `Residue`
- **CLI:** `--input proteins.fasta --enzyme Trypsin --missed-cleavages 1 --output peptide_properties.tsv`
- **Inputs:** FASTA file, enzyme, missed cleavages
- **Outputs:** TSV (peptide, protein, mass, length, charge_at_ph7, gravy, is_proteotypic, detectability_score)
- **Domain:** proteomics
- **Why needed:** Critical for targeted assay design and understanding protein coverage gaps. No TOPP tool.
- **Score:** P=7, U=8, V=7 → **7**

### 24. spectral_counting_quantifier
- **Description:** Calculate semi-quantitative protein abundances from spectral counts using emPAI, NSAF, or SIN methods.
- **pyopenms classes:** `FASTAFile`, `PeptideIdentification`, `ProteinIdentification`, `ProteaseDigestion`
- **CLI:** `--input identifications.idXML --fasta database.fasta --method nsaf --output abundances.tsv`
- **Inputs:** Identification results (idXML), FASTA database, method choice
- **Outputs:** TSV (accession, spectral_count, unique_peptides, empai/nsaf/sin_score)
- **Domain:** proteomics
- **Why needed:** Simplest label-free quantification approach. Crux toolkit has `spectral-counts` but no Python CLI with flexible input. No TOPP tool.
- **Score:** P=7, U=9, V=7 → **8**

### 25. ptm_site_localization_scorer
- **Description:** Score phosphorylation (or other PTM) site localization confidence using an Ascore-like approach comparing ion coverage.
- **pyopenms classes:** `TheoreticalSpectrumGenerator`, `SpectrumAlignment`, `AASequence`, `MSSpectrum`
- **CLI:** `--spectrum scan.mzML --scan-index 500 --peptide "PEPS(Phospho)TIDEK" --output scores.tsv`
- **Inputs:** mzML spectrum + scan, modified peptide with ambiguous site(s)
- **Outputs:** TSV (site_position, residue, ascore, probability, best_localization)
- **Domain:** proteomics
- **Why needed:** Phosphoproteomics requires localization scores. pyAscore exists in Python but lacks pyopenms integration. No TOPP tool for standalone scoring.
- **Score:** P=8, U=9, V=8 → **8**

### 26. transition_list_generator
- **Description:** Generate SRM/MRM/PRM transition lists from peptide sequences with configurable ion types, charge states, and product ion selection.
- **pyopenms classes:** `TargetedExperiment`, `ReactionMonitoringTransition`, `TheoreticalSpectrumGenerator`, `AASequence`, `TraMLFile`
- **CLI:** `--peptides PEPTIDEK,ANOTHERPEPTIDE --charge 2,3 --product-ions y3-y8 --output transitions.tsv`
- **Inputs:** Peptide list or FASTA, charge states, product ion range
- **Outputs:** TraML or Skyline-compatible CSV (precursor_mz, product_mz, ion_type, peptide, protein)
- **Domain:** proteomics
- **Why needed:** MRMMapper exists but maps existing chromatograms; this GENERATES new transitions. Essential for targeted proteomics method setup.
- **Score:** P=7, U=8, V=8 → **8**

### 27. irt_calculator
- **Description:** Convert observed retention times to indexed retention time (iRT) values using reference peptides (Biognosys iRT kit or custom standards).
- **pyopenms classes:** `PeptideIdentification`, `TransformationDescription`
- **CLI:** `--input identifications.idXML --reference-peptides irt_standards.tsv --output irt_converted.tsv`
- **Inputs:** Identification results with RT, reference peptide list (sequence, expected_iRT)
- **Outputs:** TSV (peptide, observed_rt, irt, r_squared), linear fit parameters
- **Domain:** proteomics
- **Why needed:** iRT values enable cross-lab RT comparison. OpenSwathRTNormalizer works on chromatograms, not identifications. This operates on ID-level data. No simple TOPP equivalent.
- **Score:** P=7, U=8, V=7 → **7**

### 28. rt_prediction_additive
- **Description:** Predict peptide retention time using additive amino acid hydrophobicity models (Krokhin, Guo '86, custom coefficients).
- **pyopenms classes:** `AASequence`, `Residue`
- **CLI:** `--sequence PEPTIDEK --model krokhin` or `--input peptides.tsv --output rt_predictions.tsv`
- **Inputs:** Peptide sequence(s), model choice
- **Outputs:** TSV (sequence, predicted_rt, model_used)
- **Domain:** proteomics
- **Why needed:** pyteomics `achrom` module provides this (151 stars). A pyopenms equivalent fills a gap. No TOPP tool for additive RT prediction.
- **Score:** P=6, U=8, V=6 → **7**

### 29. isoelectric_point_calculator
- **Description:** Calculate isoelectric point (pI) for peptides/proteins using Henderson-Hasselbalch with configurable pK sets (Lehninger, Sillero, Dawson).
- **pyopenms classes:** `AASequence`, `Residue`
- **CLI:** `--sequence ACDEFGHIKLMNPQRSTVWY --pk-set lehninger` or `--fasta proteins.fasta --output pi_values.tsv`
- **Inputs:** Sequence(s) or FASTA, pK set choice
- **Outputs:** TSV (sequence, pI, charge_at_pH7, mw)
- **Domain:** proteomics
- **Why needed:** pyteomics electrochem module does this. No pyopenms/TOPP equivalent. Important for 2D gel analysis, IEF optimization, and peptide fractionation.
- **Score:** P=7, U=8, V=7 → **7**

### 30. amino_acid_composition_analyzer
- **Description:** Analyze amino acid frequency and composition statistics for proteins in a FASTA file or peptide list.
- **pyopenms classes:** `FASTAFile`, `AASequence`
- **CLI:** `--input proteins.fasta --output aa_composition.tsv --per-protein`
- **Inputs:** FASTA file or peptide list
- **Outputs:** TSV (accession/peptide, A_count, C_count, ..., Y_count, total_residues, mw)
- **Domain:** proteomics
- **Why needed:** Basic bioinformatics utility needed for bias analysis, labeling efficiency estimation (e.g., counting K/R for TMT). No TOPP tool.
- **Score:** P=5, U=7, V=6 → **6**

### 31. charge_state_predictor
- **Description:** Predict expected charge state distribution for peptides based on the number of basic residues (K, R, H) and N-terminal charge.
- **pyopenms classes:** `AASequence`, `Residue`
- **CLI:** `--sequence PEPTIDEK --ph 2.0 --output charges.json` or `--input peptides.tsv --output charge_predictions.tsv`
- **Inputs:** Peptide sequence(s), pH
- **Outputs:** TSV/JSON (sequence, min_charge, max_charge, most_likely_charge, basic_residues)
- **Domain:** proteomics
- **Why needed:** Useful for method development and inclusion list creation. No TOPP tool.
- **Score:** P=5, U=8, V=6 → **6**

---

## Category 3: FASTA Database Tools (8 tools)

### 32. fasta_subset_extractor
- **Description:** Extract proteins from FASTA by accession list, header keywords, taxonomy patterns, or sequence length range.
- **pyopenms classes:** `FASTAFile`
- **CLI:** `--input uniprot.fasta --accessions accession_list.txt --output subset.fasta` or `--keyword "Homo sapiens" --min-length 50`
- **Inputs:** FASTA file, filter criteria (accession list, keyword, length range)
- **Outputs:** Filtered FASTA file
- **Domain:** proteomics
- **Why needed:** One of the most common BioStars questions. Galaxy tutorials devote entire sections to this. fasta_utilities (GitHub) addresses it. No TOPP tool.
- **Score:** P=9, U=8, V=9 → **9**

### 33. fasta_statistics_reporter
- **Description:** Report comprehensive statistics for a FASTA database: protein count, length distribution, amino acid frequency, tryptic peptide count, duplicate detection.
- **pyopenms classes:** `FASTAFile`, `ProteaseDigestion`, `AASequence`
- **CLI:** `--input database.fasta --enzyme Trypsin --output stats.json`
- **Inputs:** FASTA file, optional enzyme for peptide counting
- **Outputs:** JSON (n_proteins, avg_length, min_length, max_length, aa_frequencies, n_tryptic_peptides, n_duplicates)
- **Domain:** proteomics
- **Why needed:** Corrupt or duplicate-laden databases cause silent search errors. No TOPP tool reports these statistics comprehensively.
- **Score:** P=7, U=9, V=8 → **8**

### 34. contaminant_database_merger
- **Description:** Download/append cRAP contaminant sequences to a FASTA, add configurable prefix (e.g., CONT_), remove duplicate sequences, standardize headers.
- **pyopenms classes:** `FASTAFile`
- **CLI:** `--input target.fasta --add-crap --prefix CONT_ --remove-duplicates --output merged.fasta`
- **Inputs:** Target FASTA, contaminant addition options
- **Outputs:** Merged FASTA with contaminants prefixed
- **Domain:** proteomics
- **Why needed:** Standard database preparation step. Galaxy tutorials cover it. cRAP is from thegpm.org. No TOPP tool for this specific merge + dedup + prefix workflow.
- **Score:** P=8, U=8, V=9 → **8**

### 35. fasta_cleaner
- **Description:** Clean and validate FASTA files: remove duplicate sequences, fix header formatting, filter by length, remove stop codons and non-standard amino acids.
- **pyopenms classes:** `FASTAFile`
- **CLI:** `--input messy.fasta --remove-duplicates --fix-headers --remove-stop-codons --min-length 6 --output clean.fasta`
- **Inputs:** FASTA file, cleaning options
- **Outputs:** Cleaned FASTA file, report of changes made
- **Domain:** proteomics
- **Why needed:** Downloaded databases often have issues. profasta (GitHub) partially addresses this. No TOPP tool.
- **Score:** P=6, U=8, V=7 → **7**

### 36. fasta_merger
- **Description:** Merge multiple FASTA files with duplicate removal and header conflict resolution. Useful for combining databases from different sources.
- **pyopenms classes:** `FASTAFile`
- **CLI:** `--inputs db1.fasta db2.fasta db3.fasta --remove-duplicates --output merged.fasta`
- **Inputs:** Multiple FASTA files
- **Outputs:** Merged FASTA with duplicate handling report
- **Domain:** proteomics
- **Why needed:** Common when combining UniProt + custom sequences + contaminants. No TOPP tool for multi-file merge with dedup.
- **Score:** P=6, U=7, V=7 → **7**

### 37. fasta_decoy_validator
- **Description:** Check if a FASTA database already contains decoys, report target/decoy ratio, validate decoy prefix consistency.
- **pyopenms classes:** `FASTAFile`
- **CLI:** `--input database.fasta --decoy-prefix DECOY_ --output validation.json`
- **Inputs:** FASTA file, expected decoy prefix
- **Outputs:** JSON (has_decoys, n_target, n_decoy, ratio, prefix_consistent, issues)
- **Domain:** proteomics
- **Why needed:** Running a search with a database that already has decoys + adding more causes inflated FDR. Common mistake. No TOPP tool validates this.
- **Score:** P=6, U=9, V=7 → **7**

### 38. fasta_in_silico_digest_stats
- **Description:** Digest a FASTA database and report peptide-level statistics: total unique peptides, mass distribution, length distribution, peptides per protein.
- **pyopenms classes:** `FASTAFile`, `ProteaseDigestion`, `AASequence`
- **CLI:** `--input database.fasta --enzyme Trypsin --missed-cleavages 2 --min-length 7 --max-length 50 --output digest_stats.tsv`
- **Inputs:** FASTA, enzyme, digestion parameters
- **Outputs:** TSV (peptide stats), JSON (summary with distributions)
- **Domain:** proteomics
- **Why needed:** Understanding the search space size and peptide mass distribution is important for search parameter optimization. No TOPP tool outputs this.
- **Score:** P=6, U=8, V=7 → **7**

### 39. fasta_taxonomy_splitter
- **Description:** Split a multi-organism FASTA file by taxonomy (parsed from headers) into separate per-species files.
- **pyopenms classes:** `FASTAFile`
- **CLI:** `--input combined.fasta --header-pattern "OS=([^=]+) OX=" --output-dir split/`
- **Inputs:** FASTA file, header taxonomy pattern
- **Outputs:** Multiple FASTA files (one per species)
- **Domain:** proteomics
- **Why needed:** Metaproteomics and multi-species databases need splitting for analysis. No TOPP tool.
- **Score:** P=5, U=8, V=6 → **6**

---

## Category 4: File Conversion & Export (10 tools)

### 40. mzml_to_mgf_converter
- **Description:** Convert MS2 spectra from mzML to MGF format, preserving precursor m/z, charge, RT. Simpler and more configurable than FileConverter for this specific task (title formatting, missing charge handling).
- **pyopenms classes:** `MSExperiment`, `MzMLFile`, `MSSpectrum`
- **CLI:** `--input run.mzML --ms-level 2 --title-format "scan={scan}" --output spectra.mgf`
- **Inputs:** mzML file, MS level, title format template
- **Outputs:** MGF file
- **Domain:** proteomics / metabolomics
- **Why needed:** Most common format conversion need. FileConverter exists but this adds value: configurable titles, charge imputation for missing charges, scan number formatting for GNPS compatibility.
- **Score:** P=9, U=7, V=9 → **8**

### 41. ms_data_to_csv_exporter
- **Description:** Export MS data (spectra peaks, feature maps, consensus maps) to flat CSV/TSV. More flexible than TextExporter with column selection and filtering.
- **pyopenms classes:** `MSExperiment`, `FeatureMap`, `ConsensusMap`, `FeatureXMLFile`, `ConsensusXMLFile`
- **CLI:** `--input features.featureXML --columns mz,rt,intensity,charge,quality --output features.tsv` or `--input run.mzML --type peaks --output peaks.tsv`
- **Inputs:** mzML, featureXML, or consensusXML
- **Outputs:** TSV with selected columns
- **Domain:** proteomics / metabolomics
- **Why needed:** TextExporter is limited in column selection and doesn't handle all formats well. pandas df export from pyopenms is documented but not CLI-wrapped.
- **Score:** P=7, U=7, V=8 → **7**

### 42. mztab_summarizer
- **Description:** Parse mzTab files and extract summary statistics: protein/peptide/PSM counts, FDR stats, quantification overview, modification frequencies.
- **pyopenms classes:** `MzTabFile`, `MzTab`
- **CLI:** `--input results.mzTab --output summary.tsv`
- **Inputs:** mzTab file
- **Outputs:** TSV summary (n_proteins, n_peptides, n_psms, fdr_stats, top_modifications, quant_summary)
- **Domain:** proteomics / metabolomics
- **Why needed:** mzTab is the HUPO-PSI standard exchange format for ProteomeXchange/PRIDE but hard to read manually. No TOPP tool summarizes it.
- **Score:** P=7, U=9, V=7 → **8**

### 43. consensus_map_to_matrix
- **Description:** Convert a consensusXML file to a flat quantification matrix (rows=features, columns=samples) suitable for statistical analysis.
- **pyopenms classes:** `ConsensusMap`, `ConsensusXMLFile`
- **CLI:** `--input consensus.consensusXML --output quant_matrix.tsv --include-metadata`
- **Inputs:** ConsensusXML file
- **Outputs:** TSV matrix (feature_id, mz, rt, sample1_intensity, sample2_intensity, ...)
- **Domain:** proteomics / metabolomics
- **Why needed:** Researchers need flat matrices for R/Python statistical tools. TextExporter output is not matrix-shaped. This bridges the gap.
- **Score:** P=7, U=8, V=8 → **8**

### 44. idxml_to_tsv_exporter
- **Description:** Export idXML identification results to a flat TSV with configurable columns (peptide, protein, score, mz, rt, charge, modifications, q-value).
- **pyopenms classes:** `IdXMLFile`, `PeptideIdentification`, `ProteinIdentification`
- **CLI:** `--input results.idXML --columns peptide,protein,score,mz,rt,charge,modifications --output results.tsv`
- **Inputs:** idXML file, column selection
- **Outputs:** Flat TSV
- **Domain:** proteomics
- **Why needed:** idXML is XML and hard to work with directly. TextExporter has limited configurability. This provides a clean export.
- **Score:** P=7, U=7, V=8 → **7**

### 45. featurexml_merger
- **Description:** Merge multiple featureXML files from different tools/runs into a single featureXML, annotating source.
- **pyopenms classes:** `FeatureMap`, `FeatureXMLFile`
- **CLI:** `--inputs run1.featureXML run2.featureXML --output merged.featureXML`
- **Inputs:** Multiple featureXML files
- **Outputs:** Merged featureXML with source annotation
- **Domain:** proteomics / metabolomics
- **Why needed:** Useful when combining features from different algorithms. No TOPP tool for featureXML merging (FeatureLinker links, doesn't merge).
- **Score:** P=5, U=8, V=6 → **6**

### 46. mgf_to_mzml_converter
- **Description:** Convert MGF files back to mzML format with proper metadata (instrument, source file, etc.). Reverse of tool #40.
- **pyopenms classes:** `MSExperiment`, `MzMLFile`, `MSSpectrum`
- **CLI:** `--input spectra.mgf --output spectra.mzML`
- **Inputs:** MGF file
- **Outputs:** mzML file
- **Domain:** proteomics / metabolomics
- **Why needed:** Some tools output MGF but downstream tools need mzML. FileConverter handles this but has issues with metadata preservation.
- **Score:** P=6, U=6, V=7 → **6**

### 47. sirius_exporter
- **Description:** Export feature and MS2 data into SIRIUS-compatible .ms format for molecular formula prediction and compound class annotation.
- **pyopenms classes:** `FeatureMap`, `MSExperiment`, `FeatureXMLFile`
- **CLI:** `--features features.featureXML --mzml data.mzML --output sirius_input.ms`
- **Inputs:** featureXML + mzML
- **Outputs:** SIRIUS .ms format file
- **Domain:** metabolomics
- **Why needed:** SIRIUS is the most popular metabolomics annotation tool. UmetaFlow (published, J Cheminformatics) includes this step. No standalone CLI.
- **Score:** P=7, U=8, V=7 → **7**

### 48. spectral_library_builder
- **Description:** Build a consensus spectral library from DDA identifications + spectra. Merge replicate PSMs for same peptide, filter by FDR, output MSP/TraML.
- **pyopenms classes:** `MSExperiment`, `MzMLFile`, `IdXMLFile`, `PeptideIdentification`, `TraMLFile`
- **CLI:** `--spectra run.mzML --identifications results.idXML --fdr 0.01 --output library.msp`
- **Inputs:** mzML + idXML, FDR threshold
- **Outputs:** Spectral library (MSP, TraML, or TSV)
- **Domain:** proteomics
- **Why needed:** DIA workflows are growing rapidly. Building libraries from DDA data requires merging + filtering. No single TOPP tool does end-to-end library building.
- **Score:** P=8, U=9, V=8 → **8**

### 49. psm_feature_extractor
- **Description:** Extract rescoring features from PSMs: precursor mass error, isotope pattern fit, fragment ion coverage, matched/total intensity ratio, delta score.
- **pyopenms classes:** `MSExperiment`, `PeptideIdentification`, `TheoreticalSpectrumGenerator`, `SpectrumAlignment`
- **CLI:** `--mzml run.mzML --identifications results.idXML --output features.tsv`
- **Inputs:** mzML + idXML
- **Outputs:** TSV (scan, peptide, mass_error_ppm, fragment_coverage, matched_intensity_ratio, delta_score, ...)
- **Domain:** proteomics
- **Why needed:** MS2Rescore (63 stars) and Mokapot (50 stars) need these features. No standalone feature extraction CLI. No TOPP tool.
- **Score:** P=7, U=9, V=7 → **8**

---

## Category 5: Quality Control & Metrics (8 tools)

### 50. lc_ms_qc_reporter
- **Description:** Generate comprehensive QC report from mzML: MS1/MS2 counts, TIC stability, precursor charge distribution, mass accuracy distribution, injection time stats, MS2 trigger rate.
- **pyopenms classes:** `MSExperiment`, `MzMLFile`, `MSSpectrum`
- **CLI:** `--input run.mzML --output qc_report.json --format json`
- **Inputs:** mzML file, optional identifications (idXML)
- **Outputs:** JSON/TSV QC report with all metrics
- **Domain:** proteomics / metabolomics
- **Why needed:** QCCalculator TOPP tool exists but is limited and outputs mzQC only. This provides a more comprehensive, human-readable report.
- **Score:** P=8, U=7, V=9 → **8**

### 51. mzqc_generator
- **Description:** Generate mzQC-format (HUPO-PSI standard) quality control files from mzML data. Implements standard QC metric vocabulary.
- **pyopenms classes:** `MSExperiment`, `MzMLFile`
- **CLI:** `--input run.mzML --identifications results.idXML --output qc_report.mzQC`
- **Inputs:** mzML, optional idXML
- **Outputs:** mzQC JSON file (HUPO-PSI standard)
- **Domain:** proteomics / metabolomics
- **Why needed:** mzQC (32 stars, new HUPO-PSI standard) is growing but generation tools are scarce. QCCalculator has limited metrics. pymzqc library helps but no CLI.
- **Score:** P=7, U=8, V=8 → **8**

### 52. identification_qc_reporter
- **Description:** Report identification-level QC metrics from idXML: PSM/peptide/protein counts at FDR thresholds, score distributions, missed cleavage rates, modification frequencies, mass error distribution.
- **pyopenms classes:** `IdXMLFile`, `PeptideIdentification`, `ProteinIdentification`
- **CLI:** `--input results.idXML --fdr-thresholds 0.01,0.05 --output id_qc.json`
- **Inputs:** idXML file, FDR thresholds
- **Outputs:** JSON (n_psms, n_peptides, n_proteins per FDR, score_distribution, mass_error_stats, mc_distribution, mod_frequencies)
- **Domain:** proteomics
- **Why needed:** Identification statistics are manually computed. No TOPP tool provides this comprehensive summary.
- **Score:** P=7, U=9, V=8 → **8**

### 53. run_comparison_reporter
- **Description:** Compare two or more mzML files side-by-side: TIC overlap, shared precursors, RT shift, intensity correlation.
- **pyopenms classes:** `MSExperiment`, `MzMLFile`
- **CLI:** `--inputs run1.mzML run2.mzML --output comparison.json`
- **Inputs:** Two or more mzML files
- **Outputs:** JSON (n_spectra_diff, rt_range_diff, tic_correlation, shared_precursors, unique_precursors_per_run)
- **Domain:** proteomics / metabolomics
- **Why needed:** Comparing technical replicates or runs before/after instrument maintenance. No TOPP tool.
- **Score:** P=6, U=9, V=7 → **7**

### 54. mass_error_distribution_analyzer
- **Description:** Compute and report precursor and fragment mass error distributions from identification results. Essential for calibration decisions.
- **pyopenms classes:** `IdXMLFile`, `PeptideIdentification`, `MSExperiment`
- **CLI:** `--input results.idXML --mzml run.mzML --output mass_errors.tsv`
- **Inputs:** idXML + mzML
- **Outputs:** TSV (scan, peptide, precursor_error_ppm, fragment_errors), summary stats (median, mean, std, percentiles)
- **Domain:** proteomics
- **Why needed:** Understanding mass error distribution informs search tolerance settings and calibration needs. No TOPP tool outputs this analysis.
- **Score:** P=7, U=8, V=8 → **8**

### 55. acquisition_rate_analyzer
- **Description:** Analyze MS1/MS2 acquisition rates over time: scans per second, cycle time, duty cycle, idle time estimation.
- **pyopenms classes:** `MSExperiment`, `MzMLFile`, `MSSpectrum`
- **CLI:** `--input run.mzML --output acquisition_rate.tsv`
- **Inputs:** mzML file
- **Outputs:** TSV (time_window, ms1_rate, ms2_rate, cycle_time, duty_cycle)
- **Domain:** proteomics
- **Why needed:** Important for method optimization and comparing DDA vs DIA efficiency. No TOPP tool.
- **Score:** P=6, U=9, V=7 → **7**

### 56. precursor_isolation_purity
- **Description:** Estimate precursor isolation purity by examining the MS1 spectrum around each selected precursor. Reports co-isolation interference.
- **pyopenms classes:** `MSExperiment`, `MSSpectrum`, `Precursor`
- **CLI:** `--input run.mzML --output purity.tsv`
- **Inputs:** mzML file
- **Outputs:** TSV (ms2_scan, precursor_mz, isolation_width, purity_fraction, n_interfering_peaks)
- **Domain:** proteomics
- **Why needed:** Co-isolation is a major source of quantification error in DDA/TMT. Chimeric spectra reduce identification rates. No TOPP tool computes this metric.
- **Score:** P=7, U=9, V=7 → **8**

### 57. injection_time_analyzer
- **Description:** Extract and analyze injection time (ion accumulation time) values from mzML metadata across all scans.
- **pyopenms classes:** `MSExperiment`, `MSSpectrum`
- **CLI:** `--input run.mzML --output injection_times.tsv`
- **Inputs:** mzML file
- **Outputs:** TSV (scan, ms_level, rt, injection_time_ms), summary stats
- **Domain:** proteomics
- **Why needed:** Injection time is crucial for dynamic range and sensitivity analysis in Orbitrap data. Not extracted by any TOPP tool as standalone output.
- **Score:** P=5, U=9, V=6 → **7**

---

## Category 6: Metabolomics-Specific Tools (12 tools)

### 58. adduct_calculator
- **Description:** Given a molecular formula or exact mass, compute expected m/z for all common ESI adducts in positive and/or negative mode.
- **pyopenms classes:** `EmpiricalFormula`, `Element`
- **CLI:** `--formula C6H12O6 --mode positive --output adducts.tsv` or `--mass 180.0634 --mode both`
- **Inputs:** Molecular formula or exact mass, ionization mode
- **Outputs:** TSV (adduct_type, charge, expected_mz, delta_mass)
- **Domain:** metabolomics
- **Why needed:** Fiehn Lab MS Adduct Calculator is one of the most used metabolomics web tools. MSAC paper cited 40+ times. No offline CLI. No TOPP tool.
- **Score:** P=9, U=9, V=9 → **9**

### 59. metabolite_formula_annotator
- **Description:** Annotate features with candidate molecular formulas based on accurate mass and isotope pattern matching.
- **pyopenms classes:** `EmpiricalFormula`, `CoarseIsotopePatternGenerator`, `FineIsotopePatternGenerator`
- **CLI:** `--input features.tsv --ppm 5 --elements C,H,N,O,S,P --output annotated.tsv`
- **Inputs:** Feature list (m/z, intensity, optional isotope peaks), ppm tolerance
- **Outputs:** TSV (mz, candidate_formulas, scores, isotope_fit)
- **Domain:** metabolomics
- **Why needed:** Links mass decomposition with isotope pattern scoring. No TOPP tool combines both steps.
- **Score:** P=7, U=8, V=8 → **8**

### 60. adduct_group_analyzer
- **Description:** Given a feature list, identify groups of features that likely originate from the same metabolite (different adducts of the same compound) based on mass relationships.
- **pyopenms classes:** `EmpiricalFormula`
- **CLI:** `--input features.tsv --rt-tolerance 5 --adducts "[M+H]+,[M+Na]+,[M+K]+,[M+NH4]+" --output groups.tsv`
- **Inputs:** Feature list (mz, rt, intensity), RT tolerance, adduct list
- **Outputs:** TSV (feature_id, group_id, adduct_type, neutral_mass)
- **Domain:** metabolomics
- **Why needed:** MetaboliteAdductDecharger is TOPP but works on featureXML. This works on simple TSV input, making it accessible without the full OpenMS pipeline.
- **Score:** P=7, U=7, V=8 → **7**

### 61. isotope_pattern_scorer
- **Description:** Score how well an observed isotope pattern matches the theoretical pattern for a given molecular formula. Reports fit quality metrics.
- **pyopenms classes:** `EmpiricalFormula`, `CoarseIsotopePatternGenerator`, `FineIsotopePatternGenerator`
- **CLI:** `--observed-peaks "180.063:100,181.067:6.5,182.070:0.5" --formula C6H12O6 --output fit.json`
- **Inputs:** Observed peak list (mz:intensity pairs), candidate formula
- **Outputs:** JSON (formula, chi_squared, cosine_similarity, rms_error, peak_by_peak_comparison)
- **Domain:** metabolomics
- **Why needed:** Isotope pattern scoring is crucial for formula ranking but no standalone CLI exists. Complements isotope_pattern_matcher (already in repo) with scoring.
- **Score:** P=6, U=8, V=7 → **7**

### 62. mass_difference_network_builder
- **Description:** Build a mass difference network from features: connect features whose mass difference matches known biotransformations (e.g., +15.995 oxidation, +176.032 glucuronidation).
- **pyopenms classes:** `EmpiricalFormula`
- **CLI:** `--input features.tsv --reactions biotransformations.tsv --tolerance 0.005 --output network.tsv`
- **Inputs:** Feature list (mz, rt), biotransformation mass list, tolerance
- **Outputs:** TSV (feature1, feature2, mass_diff, reaction_name, rt_diff)
- **Domain:** metabolomics
- **Why needed:** Core technique for drug metabolism and untargeted metabolomics. GNPS molecular networking uses this principle. No TOPP tool.
- **Score:** P=7, U=9, V=7 → **8**

### 63. targeted_feature_extractor
- **Description:** Extract features for a defined set of compounds (with known formulas and expected RTs) from MS1 data. Simpler than FeatureFinderMetaboIdent for quick targeted lookups.
- **pyopenms classes:** `MSExperiment`, `MzMLFile`, `EmpiricalFormula`
- **CLI:** `--input sample.mzML --targets compounds.tsv --rt-tolerance 30 --mz-ppm 5 --output quantified.tsv`
- **Inputs:** mzML file, target compound list (name, formula, expected_rt)
- **Outputs:** TSV (compound, observed_mz, observed_rt, intensity, peak_area, found)
- **Domain:** metabolomics
- **Why needed:** Quick targeted compound extraction without full feature detection pipeline. No simple TOPP equivalent for TSV-in/TSV-out workflow.
- **Score:** P=7, U=7, V=8 → **7**

### 64. blank_subtraction_tool
- **Description:** Subtract blank/control features from sample features based on m/z and RT matching. Essential metabolomics preprocessing step.
- **pyopenms classes:** `MSExperiment`, `MzMLFile`
- **CLI:** `--sample sample.mzML --blank blank.mzML --fold-change 3 --output cleaned.tsv`
- **Inputs:** Sample mzML, blank mzML, fold-change threshold
- **Outputs:** TSV of features present in sample but not in blank (or above fold-change threshold)
- **Domain:** metabolomics
- **Why needed:** Blank subtraction is one of the first steps in untargeted metabolomics. Usually done in R/Excel. No TOPP tool.
- **Score:** P=7, U=9, V=8 → **8**

### 65. retention_index_calculator
- **Description:** Calculate Kovats retention indices using alkane standard retention times. Maps observed RTs to standardized RI values.
- **pyopenms classes:** `MSExperiment` (for chromatogram access)
- **CLI:** `--input run.mzML --alkane-standards standards.tsv --output ri_converted.tsv`
- **Inputs:** mzML (or feature list), alkane standard RT/carbon number table
- **Outputs:** TSV (feature_id, observed_rt, retention_index)
- **Domain:** metabolomics
- **Why needed:** RI standardization is essential for GC-MS metabolomics compound identification. No TOPP tool.
- **Score:** P=6, U=9, V=6 → **7**

### 66. massql_query_tool
- **Description:** Query mzML data using MassQL (Mass Spec Query Language) — a SQL-like language for MS data.
- **pyopenms classes:** `MSExperiment`, `MzMLFile`
- **CLI:** `--input data.mzML --query "QUERY scaninfo(MS2DATA) WHERE MS2PROD=226.18" --output results.tsv`
- **Inputs:** mzML file, MassQL query string
- **Outputs:** TSV of matching spectra/features
- **Domain:** metabolomics / proteomics
- **Why needed:** MassQL is documented in pyopenms docs as a powerful query tool. No TOPP tool. Enables complex spectrum filtering without custom code.
- **Score:** P=6, U=9, V=6 → **7**

### 67. metabolite_class_annotator
- **Description:** Annotate features with probable metabolite classes based on mass defect ranges, adduct patterns, and RT windows characteristic of lipids, amino acids, nucleotides, etc.
- **pyopenms classes:** `EmpiricalFormula`
- **CLI:** `--input features.tsv --class-rules metabolite_classes.tsv --output annotated.tsv`
- **Inputs:** Feature list, classification rules (TSV defining mass ranges, defect ranges per class)
- **Outputs:** TSV (feature_id, mz, rt, predicted_class, confidence)
- **Domain:** metabolomics
- **Why needed:** Quick class-level annotation before running SIRIUS/CANOPUS. No TOPP tool.
- **Score:** P=6, U=9, V=6 → **7**

### 68. duplicate_feature_detector
- **Description:** Detect and flag duplicate/redundant features in a feature list based on m/z and RT proximity. Common artifact in feature detection.
- **pyopenms classes:** `FeatureMap`, `FeatureXMLFile` (optional — can work on TSV)
- **CLI:** `--input features.tsv --mz-tolerance 10ppm --rt-tolerance 5 --output deduplicated.tsv`
- **Inputs:** Feature list (mz, rt, intensity)
- **Outputs:** TSV with duplicate groups annotated, deduplicated feature list
- **Domain:** metabolomics
- **Why needed:** Feature detection often produces redundant features. Manual deduplication is error-prone. No TOPP tool.
- **Score:** P=6, U=8, V=7 → **7**

### 69. formula_mass_calculator
- **Description:** Calculate exact monoisotopic/average masses for molecular formulas, with support for adducts and charge states. Batch mode for formula lists.
- **pyopenms classes:** `EmpiricalFormula`, `Element`
- **CLI:** `--formula C6H12O6 --adduct "[M+H]+" --output mass.json` or `--formula-list formulas.tsv --output masses.tsv`
- **Inputs:** Molecular formula(s), optional adduct, charge
- **Outputs:** JSON/TSV (formula, mono_mass, avg_mass, mz, elements)
- **Domain:** metabolomics
- **Why needed:** Quick scriptable mass calculations for metabolomics. No TOPP tool for formula-to-mass with adduct support.
- **Score:** P=7, U=7, V=8 → **7**

---

## Category 7: RNA / Oligonucleotide Tools (3 tools)

### 70. rna_mass_calculator
- **Description:** Calculate monoisotopic/average mass, molecular formula, and isotope distribution for RNA/oligonucleotide sequences with modifications.
- **pyopenms classes:** `NASequence`, `RibonucleotideDB`, `EmpiricalFormula`
- **CLI:** `--sequence AAUGCAAUGG --charge 3 --modifications "[m1A]" --output rna_mass.json`
- **Inputs:** RNA sequence, charge, optional modifications
- **Outputs:** JSON (sequence, mono_mass, avg_mass, mz, formula, isotope_distribution)
- **Domain:** proteomics (nucleic acids)
- **Why needed:** Growing field of RNA therapeutics and oligonucleotide MS. No TOPP tool for RNA mass calculation.
- **Score:** P=6, U=9, V=6 → **7**

### 71. rna_digest
- **Description:** Perform in silico digestion of RNA sequences using various RNases (RNase T1, RNase U2, etc.).
- **pyopenms classes:** `RNaseDigestion`, `RNaseDB`, `NASequence`
- **CLI:** `--sequence AAUGCAAUGG --enzyme RNase_T1 --missed-cleavages 1 --output fragments.tsv`
- **Inputs:** RNA sequence, enzyme, missed cleavage limit
- **Outputs:** TSV (fragment_sequence, start, end, mono_mass, charge_states)
- **Domain:** proteomics (nucleic acids)
- **Why needed:** RNA MS is growing. No TOPP tool for RNA digestion.
- **Score:** P=5, U=9, V=5 → **6**

### 72. rna_fragment_spectrum_generator
- **Description:** Generate theoretical fragment spectra for RNA/oligonucleotide sequences (c/y/w/a-B ion series).
- **pyopenms classes:** `NASequence`, `TheoreticalSpectrumGenerator` (RNA mode), `MSSpectrum`
- **CLI:** `--sequence AAUGC --charge 2 --output rna_fragments.tsv`
- **Inputs:** RNA sequence, charge
- **Outputs:** TSV (ion_type, position, mz, annotation)
- **Domain:** proteomics (nucleic acids)
- **Why needed:** Complementary to peptide theoretical spectrum generator for the RNA world. No TOPP tool.
- **Score:** P=5, U=9, V=5 → **6**

---

## Category 8: Statistical & Quantification Tools (7 tools)

### 73. missing_value_imputation
- **Description:** Impute missing values in protein/peptide quantification matrices using MinDet, MinProb, KNN, or left-censored methods.
- **pyopenms classes:** Minimal — reads ConsensusXML optionally via `ConsensusXMLFile`; primarily numpy/scipy
- **CLI:** `--input quant_matrix.tsv --method knn --output imputed.tsv`
- **Inputs:** Quantification matrix (TSV), method choice
- **Outputs:** Imputed TSV matrix
- **Domain:** proteomics / metabolomics
- **Why needed:** Missing values are pervasive in LFQ (30-50% missing). Most tools are in R (DEP, MSnbase). No Python CLI. No TOPP tool.
- **Score:** P=7, U=9, V=8 → **8**

### 74. quantification_normalizer
- **Description:** Normalize a quantification matrix using median centering, quantile normalization, variance stabilization, or total intensity normalization.
- **pyopenms classes:** Minimal — can read ConsensusXML; primarily numpy/scipy
- **CLI:** `--input quant_matrix.tsv --method median --output normalized.tsv`
- **Inputs:** Quantification matrix (TSV), method
- **Outputs:** Normalized TSV matrix
- **Domain:** proteomics / metabolomics
- **Why needed:** ConsensusMapNormalizer works on consensusXML only. This works on flat TSV matrices from any source. No TOPP equivalent.
- **Score:** P=7, U=7, V=8 → **7**

### 75. differential_expression_tester
- **Description:** Run simple differential expression analysis (t-test, Welch's t-test, Mann-Whitney) on quantification matrices. Report fold changes, p-values, adjusted p-values.
- **pyopenms classes:** Minimal — reads matrices
- **CLI:** `--input quant_matrix.tsv --design experimental_design.tsv --test ttest --correction bh --output de_results.tsv`
- **Inputs:** Quantification matrix, experimental design (sample-to-group mapping)
- **Outputs:** TSV (protein/feature, fold_change, log2fc, pvalue, adj_pvalue)
- **Domain:** proteomics / metabolomics
- **Why needed:** AlphaPeptStats (86 stars) does this but is a heavy install. Simple t-test + multiple testing correction as a CLI fills a real gap.
- **Score:** P=7, U=8, V=7 → **7**

### 76. volcano_plot_data_generator
- **Description:** Generate volcano plot data (log2 fold change vs -log10 p-value) from quantification results with significance thresholds and protein/feature labeling.
- **pyopenms classes:** Minimal
- **CLI:** `--input de_results.tsv --fc-threshold 1.0 --pvalue-threshold 0.05 --output volcano_data.tsv`
- **Inputs:** Differential expression results (TSV with fold_change, pvalue columns)
- **Outputs:** TSV (protein, log2fc, neg_log10_pvalue, significant, label)
- **Domain:** proteomics / metabolomics
- **Why needed:** bioinfokit (367 stars) provides volcano plots. This generates the data for any plotting tool. No TOPP tool.
- **Score:** P=6, U=7, V=6 → **6**

### 77. sample_correlation_calculator
- **Description:** Calculate pairwise Pearson/Spearman correlations between samples in a quantification matrix. Useful for QC and outlier detection.
- **pyopenms classes:** Minimal
- **CLI:** `--input quant_matrix.tsv --method pearson --output correlations.tsv`
- **Inputs:** Quantification matrix
- **Outputs:** TSV correlation matrix
- **Domain:** proteomics / metabolomics
- **Why needed:** Basic QC metric for multi-sample experiments. Usually done in R. No TOPP tool.
- **Score:** P=5, U=7, V=6 → **6**

### 78. coefficient_of_variation_calculator
- **Description:** Calculate CV% (coefficient of variation) for features across replicates. Standard QC metric for quantification reproducibility.
- **pyopenms classes:** Minimal
- **CLI:** `--input quant_matrix.tsv --groups replicate_groups.tsv --output cv_report.tsv`
- **Inputs:** Quantification matrix, replicate grouping
- **Outputs:** TSV (feature, mean, sd, cv_percent), summary (median_cv, features_below_20pct_cv)
- **Domain:** proteomics / metabolomics
- **Why needed:** CV calculation is a universal QC metric but usually done in spreadsheets. No TOPP tool.
- **Score:** P=5, U=8, V=7 → **7**

### 79. intensity_distribution_reporter
- **Description:** Report intensity distribution statistics per sample in a quantification matrix: median, IQR, dynamic range, number of quantified features.
- **pyopenms classes:** Minimal
- **CLI:** `--input quant_matrix.tsv --output intensity_stats.tsv`
- **Inputs:** Quantification matrix
- **Outputs:** TSV (sample, n_quantified, median_intensity, iqr, dynamic_range_log10)
- **Domain:** proteomics / metabolomics
- **Why needed:** Quick per-sample intensity QC. No TOPP tool.
- **Score:** P=5, U=8, V=6 → **6**

---

## Category 9: Data Integration & Interoperability (7 tools)

### 80. search_result_merger
- **Description:** Merge identification results from multiple search engines (idXML files) into a combined result with consensus scoring.
- **pyopenms classes:** `IdXMLFile`, `PeptideIdentification`, `ProteinIdentification`
- **CLI:** `--inputs comet.idXML msgf.idXML xtandem.idXML --output consensus.idXML --method intersection`
- **Inputs:** Multiple idXML files, merge method (union/intersection)
- **Outputs:** Merged idXML with source engine annotations
- **Domain:** proteomics
- **Why needed:** Multi-engine searching improves identifications. Ursgal (45 stars) does this. No simple pyopenms CLI.
- **Score:** P=6, U=8, V=7 → **7**

### 81. peptide_to_protein_mapper
- **Description:** Map a list of peptide sequences to their parent proteins in a FASTA database. Report all matching proteins per peptide.
- **pyopenms classes:** `FASTAFile`, `AASequence`
- **CLI:** `--peptides peptide_list.tsv --fasta database.fasta --output mapped.tsv`
- **Inputs:** Peptide list, FASTA database
- **Outputs:** TSV (peptide, protein_accessions, positions, is_unique)
- **Domain:** proteomics
- **Why needed:** Reverse mapping peptides to proteins without running a full search. Useful for validation and targeted proteomics. No TOPP tool.
- **Score:** P=6, U=8, V=7 → **7**

### 82. inclusion_list_generator
- **Description:** Generate instrument inclusion/exclusion lists from identification results or feature lists for targeted re-analysis.
- **pyopenms classes:** `PeptideIdentification`, `AASequence`
- **CLI:** `--input results.idXML --format thermo --min-score 0.95 --charge 2,3 --output inclusion.csv`
- **Inputs:** idXML or feature list, instrument format, filters
- **Outputs:** Instrument-specific inclusion list CSV (precursor_mz, charge, rt_start, rt_end)
- **Domain:** proteomics
- **Why needed:** Researchers manually create inclusion lists from results. Tedious and error-prone. No TOPP tool.
- **Score:** P=7, U=9, V=7 → **8**

### 83. maxquant_result_converter
- **Description:** Convert MaxQuant output files (evidence.txt, proteinGroups.txt, msms.txt) to pyopenms-compatible formats (idXML, featureXML).
- **pyopenms classes:** `IdXMLFile`, `PeptideIdentification`, `FeatureMap`
- **CLI:** `--input evidence.txt --fasta database.fasta --output results.idXML` or `--input proteinGroups.txt --output proteins.tsv`
- **Inputs:** MaxQuant output files
- **Outputs:** idXML, featureXML, or standardized TSV
- **Domain:** proteomics
- **Why needed:** MaxQuant is the most used proteomics search engine but outputs proprietary formats. No TOPP tool for MaxQuant import.
- **Score:** P=8, U=9, V=7 → **8**

### 84. diann_result_converter
- **Description:** Convert DIA-NN report files to pyopenms-compatible formats or standardized TSV for downstream analysis.
- **pyopenms classes:** `IdXMLFile`, `PeptideIdentification`
- **CLI:** `--input report.tsv --output results.idXML` or `--input report.tsv --format mztab --output results.mztab`
- **Inputs:** DIA-NN report.tsv
- **Outputs:** idXML, mzTab, or standardized TSV
- **Domain:** proteomics
- **Why needed:** DIA-NN (431 stars) is the top DIA engine. Converting its output to standard formats enables interoperability. No TOPP tool.
- **Score:** P=8, U=9, V=7 → **8**

### 85. fragpipe_result_converter
- **Description:** Convert FragPipe/MSFragger output files (psm.tsv, combined_protein.tsv) to standard formats (idXML, mzTab).
- **pyopenms classes:** `IdXMLFile`, `PeptideIdentification`
- **CLI:** `--input psm.tsv --output results.idXML`
- **Inputs:** FragPipe output files
- **Outputs:** idXML, mzTab, or standardized TSV
- **Domain:** proteomics
- **Why needed:** FragPipe/MSFragger is widely used. Format conversion to standard outputs enables downstream tool interoperability. No TOPP tool.
- **Score:** P=7, U=9, V=7 → **8**

### 86. experimental_design_generator
- **Description:** Generate OpenMS experimental design TSV from file paths and sample annotations. Required by many OpenMS workflows but tedious to create manually.
- **pyopenms classes:** `ExperimentalDesign`
- **CLI:** `--mzml-dir runs/ --conditions condition_map.tsv --output experimental_design.tsv`
- **Inputs:** Directory of mzML files, condition mapping
- **Outputs:** OpenMS experimental design TSV
- **Domain:** proteomics / metabolomics
- **Why needed:** OpenMS workflows require this file but creating it is manual and error-prone. No generator tool.
- **Score:** P=6, U=9, V=7 → **7**

---

## Category 10: Specialized & Emerging Tools (14 tools)

### 87. crosslink_mass_calculator
- **Description:** Calculate expected masses for crosslinked peptide pairs using common crosslinkers (BS3, DSS, DSSO, etc.).
- **pyopenms classes:** `AASequence`, `EmpiricalFormula`
- **CLI:** `--peptide1 PEPTIDEK --peptide2 ANOTHERPEPTIDER --crosslinker DSS --charge 3 --output crosslink_masses.tsv`
- **Inputs:** Two peptide sequences, crosslinker name, charge
- **Outputs:** TSV (crosslinked_mass, mz, crosslinker_mass, link_sites)
- **Domain:** proteomics
- **Why needed:** XL-MS is a growing structural biology technique. No simple mass calculator exists as a CLI. No TOPP tool.
- **Score:** P=6, U=9, V=6 → **7**

### 88. glycopeptide_mass_calculator
- **Description:** Calculate masses for glycosylated peptides with common N-glycan compositions (HexNAc, Hex, Fuc, NeuAc).
- **pyopenms classes:** `AASequence`, `EmpiricalFormula`
- **CLI:** `--sequence PEPTIDEK --glycan "HexNAc(2)Hex(5)Fuc(1)" --charge 3 --output glyco_masses.tsv`
- **Inputs:** Peptide sequence, glycan composition, charge
- **Outputs:** TSV (peptide_mass, glycan_mass, total_mass, mz, composition)
- **Domain:** proteomics
- **Why needed:** Glycoproteomics is rapidly growing. Manual glycan mass calculations are error-prone. No TOPP tool.
- **Score:** P=6, U=9, V=6 → **7**

### 89. immunopeptide_filter
- **Description:** Filter identification results for MHC-I or MHC-II peptides based on length range, binding motif, and amino acid composition.
- **pyopenms classes:** `PeptideIdentification`, `AASequence`, `IdXMLFile`
- **CLI:** `--input results.idXML --class-i --length-range 8-11 --allele HLA-A02:01 --output immunopeptides.tsv`
- **Inputs:** idXML, MHC class, length range, optional allele motif
- **Outputs:** TSV (peptide, protein, length, n_terminal_aa, c_terminal_aa, motif_match)
- **Domain:** proteomics (immunopeptidomics)
- **Why needed:** MHCflurry (237 stars) predicts binding but no simple filter tool exists. Growing field with cancer immunotherapy applications.
- **Score:** P=7, U=9, V=6 → **7**

### 90. semi_tryptic_peptide_finder
- **Description:** Identify semi-tryptic and non-tryptic peptides in identification results. Important for degradome analysis and sample quality assessment.
- **pyopenms classes:** `ProteaseDigestion`, `AASequence`, `PeptideIdentification`
- **CLI:** `--input results.idXML --enzyme Trypsin --output semi_tryptic.tsv`
- **Inputs:** idXML, enzyme
- **Outputs:** TSV (peptide, protein, cleavage_type [fully_tryptic/semi_tryptic/non_tryptic], n_term_ok, c_term_ok)
- **Domain:** proteomics
- **Why needed:** Semi-tryptic peptides indicate degradation or biological cleavage. No TOPP tool classifies cleavage types per PSM.
- **Score:** P=6, U=9, V=7 → **7**

### 91. collision_energy_analyzer
- **Description:** Extract and analyze collision energy (CE) values across all MS2 spectra in a DDA run. Report CE distribution and identify CE stepping patterns.
- **pyopenms classes:** `MSExperiment`, `MSSpectrum`, `Precursor`
- **CLI:** `--input run.mzML --output ce_analysis.tsv`
- **Inputs:** mzML file
- **Outputs:** TSV (scan, precursor_mz, collision_energy, ce_type), CE distribution summary
- **Domain:** proteomics
- **Why needed:** CE optimization is critical for fragmentation quality. No TOPP tool extracts and reports CE values.
- **Score:** P=5, U=9, V=6 → **7**

### 92. ms_data_ml_exporter
- **Description:** Export MS data as feature matrices suitable for machine learning. Extract features like m/z, RT, intensity, charge, peak shape metrics.
- **pyopenms classes:** `MSExperiment`, `FeatureMap`
- **CLI:** `--input features.featureXML --features mz,rt,intensity,charge,fwhm --output ml_matrix.csv`
- **Inputs:** featureXML or mzML
- **Outputs:** CSV/numpy/parquet feature matrix
- **Domain:** proteomics / metabolomics
- **Why needed:** pyopenms docs have a dedicated ML interfacing section. No TOPP tool bridges MS data to ML formats.
- **Score:** P=6, U=8, V=6 → **7**

### 93. peptide_spectral_match_validator
- **Description:** Validate individual PSMs by re-scoring: compute fragment ion coverage, precursor mass error, b/y ion ratio, explained intensity fraction.
- **pyopenms classes:** `TheoreticalSpectrumGenerator`, `SpectrumAlignment`, `MSExperiment`, `PeptideIdentification`
- **CLI:** `--mzml run.mzML --identifications results.idXML --output validation.tsv`
- **Inputs:** mzML + idXML
- **Outputs:** TSV (scan, peptide, fragment_coverage, explained_intensity, mass_error_ppm, b_y_ratio, valid)
- **Domain:** proteomics
- **Why needed:** Standalone PSM validation without full rescoring. Useful for manual curation and suspicious hit investigation. No TOPP tool.
- **Score:** P=7, U=8, V=7 → **7**

### 94. protein_group_reporter
- **Description:** Parse protein groups from idXML and generate a clean protein-level report with group membership, razor peptides, and protein scores.
- **pyopenms classes:** `IdXMLFile`, `ProteinIdentification`
- **CLI:** `--input results.idXML --output protein_groups.tsv`
- **Inputs:** idXML with protein inference results
- **Outputs:** TSV (group_id, protein_accessions, n_peptides, n_unique_peptides, n_razor_peptides, score, coverage)
- **Domain:** proteomics
- **Why needed:** Extracting clean protein group tables from idXML is cumbersome. TextExporter output is not group-oriented. No TOPP tool.
- **Score:** P=6, U=8, V=7 → **7**

### 95. mzml_metadata_extractor
- **Description:** Extract instrument metadata from mzML files: instrument model, serial number, acquisition software, source file, contact info, processing methods.
- **pyopenms classes:** `MSExperiment`, `MzMLFile`
- **CLI:** `--input run.mzML --output metadata.json`
- **Inputs:** mzML file
- **Outputs:** JSON (instrument_model, serial_number, software, source_file, acquisition_date, contact)
- **Domain:** proteomics / metabolomics
- **Why needed:** Metadata extraction is needed for data management, LIMS integration, and ProteomeXchange submission. No TOPP tool outputs just metadata.
- **Score:** P=6, U=8, V=7 → **7**

### 96. sequence_tag_generator
- **Description:** Generate de novo sequence tags from MS2 spectra — short amino acid subsequences derived from fragment ion ladders without database search.
- **pyopenms classes:** `MSSpectrum`, `Residue`, `ResidueDB`
- **CLI:** `--input spectrum.mzML --scan-index 500 --tolerance 0.02 --min-tag-length 3 --output tags.tsv`
- **Inputs:** mzML + scan index, tolerance, minimum tag length
- **Outputs:** TSV (tag_sequence, start_mass, end_mass, score, ion_series)
- **Domain:** proteomics
- **Why needed:** Sequence tags enable error-tolerant searching and de novo sequencing validation. No TOPP tool generates tags from single spectra.
- **Score:** P=6, U=9, V=6 → **7**

### 97. precursor_recurrence_analyzer
- **Description:** Analyze how often the same precursor m/z is selected for MS2 across a run (redundant MS2 triggers). Key DDA efficiency metric.
- **pyopenms classes:** `MSExperiment`, `MSSpectrum`, `Precursor`
- **CLI:** `--input run.mzML --mz-tolerance 10ppm --rt-tolerance 30 --output recurrence.tsv`
- **Inputs:** mzML file, grouping tolerances
- **Outputs:** TSV (precursor_mz, n_selections, first_rt, last_rt, total_rt_span)
- **Domain:** proteomics
- **Why needed:** Understanding MS2 redundancy helps optimize dynamic exclusion settings. No TOPP tool.
- **Score:** P=6, U=9, V=6 → **7**

### 98. peptide_mass_fingerprint
- **Description:** Generate a peptide mass fingerprint (PMF) from in silico digestion of a protein. Compare against observed masses for simple protein identification.
- **pyopenms classes:** `ProteaseDigestion`, `AASequence`, `FASTAFile`
- **CLI:** `--fasta proteins.fasta --enzyme Trypsin --accession P12345 --output fingerprint.tsv` or `--observe masses.tsv --fasta proteins.fasta --tolerance 10ppm --output matches.tsv`
- **Inputs:** FASTA + protein accession (generation mode) or observed masses + FASTA (matching mode)
- **Outputs:** TSV (peptide, mass) or TSV (protein, matched_peptides, coverage, score)
- **Domain:** proteomics
- **Why needed:** PMF is still used in MALDI-TOF workflows. No standalone Python CLI. No TOPP tool for PMF matching.
- **Score:** P=5, U=8, V=5 → **6**

### 99. sample_complexity_estimator
- **Description:** Estimate sample complexity from MS1 data: number of distinct isotope envelopes, feature density over RT, dynamic range.
- **pyopenms classes:** `MSExperiment`, `MzMLFile`
- **CLI:** `--input run.mzML --output complexity.json`
- **Inputs:** mzML file
- **Outputs:** JSON (estimated_features, peak_density_per_second, dynamic_range_log10, max_concurrent_peptides)
- **Domain:** proteomics / metabolomics
- **Why needed:** Helps decide fractionation strategy and instrument method parameters. No TOPP tool.
- **Score:** P=5, U=9, V=5 → **6**

### 100. ms1_feature_intensity_tracker
- **Description:** Track intensity of specific MS1 features (m/z + RT) across a batch of runs. Report presence/absence and intensity trends.
- **pyopenms classes:** `MSExperiment`, `MzMLFile`
- **CLI:** `--inputs run1.mzML run2.mzML run3.mzML --features targets.tsv --ppm 10 --rt-tolerance 30 --output tracking.tsv`
- **Inputs:** Multiple mzML files, target feature list (mz, expected_rt)
- **Outputs:** TSV (feature_mz, feature_rt, run1_intensity, run2_intensity, ..., cv_percent)
- **Domain:** proteomics / metabolomics
- **Why needed:** Batch-level feature monitoring for QC and longitudinal studies. No TOPP tool.
- **Score:** P=6, U=9, V=7 → **7**

---

---

## Category 11: Niche Proteomics — Paper-Level Tools (22 tools)

*These tools address specialized workflows that researchers perform in published papers. Each references specific DOIs.*

### 101. phosphosite_class_filter
- **Description:** Classify phosphosites into Class I/II/III by localization probability. Report enrichment efficiency (phospho/total ratio) and mono/di/tri-phospho distribution.
- **Paper:** Beausoleil et al. "A probability-based approach for phosphorylation analysis" (DOI: 10.1038/nbt1240); Comparing 22 Pipelines (DOI: 10.1021/acs.jproteome.9b00679)
- **pyopenms classes:** `AASequence`, `ModificationsDB`, `ResidueModification`
- **CLI:** `--input phosphosites.tsv --class1-threshold 0.75 --output classified_sites.tsv`
- **Domain:** proteomics (phosphoproteomics)
- **Score:** P=8, U=9, V=8 → **8**

### 102. phospho_motif_analyzer
- **Description:** Extract ±7aa windows around phosphosites, compute position-specific amino acid frequencies vs. proteome background, match kinase substrate motifs (PKA: RxxS, CK2: SxxE, etc.).
- **Paper:** MMFPh (DOI: 10.1093/bioinformatics/bts195); pLogo (DOI: 10.1038/nmeth.2541)
- **pyopenms classes:** `AASequence`, `FASTAFile`, `ResidueDB`
- **CLI:** `--input phosphosites.tsv --fasta proteome.fasta --window 7 --kinase-motifs kinase_db.tsv --output motifs.tsv`
- **Domain:** proteomics (phosphoproteomics)
- **Score:** P=7, U=9, V=8 → **8**

### 103. phospho_enrichment_qc
- **Description:** From search results, compute phospho-enrichment efficiency (%), singly/multiply phosphorylated fractions, and pSer/pThr/pTyr ratios.
- **Paper:** Humphrey et al. "High-throughput phosphoproteomics" (Nat Biotechnol 2015)
- **pyopenms classes:** `AASequence`, `ModificationsDB`
- **CLI:** `--input search_results.tsv --modification Phospho --output enrichment_stats.tsv`
- **Domain:** proteomics (phosphoproteomics)
- **Score:** P=7, U=9, V=7 → **8**

### 104. proteogenomics_db_builder
- **Description:** Build custom FASTA from reference proteome + VCF variants (SAVs, indels). Generate mutant protein entries with flanking residue context.
- **Paper:** moPepGen (DOI: 10.1038/s41587-025-02701-0); Proteogenomic DB from RNA-seq (DOI: 10.1186/1471-2105-15-S7-S7)
- **pyopenms classes:** `FASTAFile`, `AASequence`
- **CLI:** `--fasta reference.fasta --vcf variants.vcf --flanking 7 --output custom_db.fasta`
- **Domain:** proteomics (proteogenomics)
- **Score:** P=7, U=9, V=7 → **8**

### 105. variant_peptide_validator
- **Description:** Validate SAV peptide identifications: check isobaric mass ambiguity with canonical peptides, count variant-site-determining fragment ions.
- **Paper:** PgxSAVy (DOI: 10.1016/j.mcpro.2023.100653)
- **pyopenms classes:** `TheoreticalSpectrumGenerator`, `SpectrumAlignment`, `AASequence`
- **CLI:** `--input variant_psms.tsv --fasta reference.fasta --spectra run.mzML --min-site-ions 3 --output validated.tsv`
- **Domain:** proteomics (proteogenomics)
- **Score:** P=6, U=9, V=7 → **7**

### 106. hdx_deuterium_uptake
- **Description:** Calculate deuterium uptake from HDX-MS data: centroid mass shift, relative fractional uptake, back-exchange correction using fully deuterated control.
- **Paper:** HDX-MS Recommendations (DOI: 10.1038/s41592-019-0459-y)
- **pyopenms classes:** `AASequence`, `EmpiricalFormula`, `MSSpectrum`
- **CLI:** `--peptides peptide_list.tsv --undeuterated ref.tsv --fully-deuterated fd.tsv --timepoints 0,10,60,300,3600 --output uptake.tsv`
- **Domain:** proteomics (structural)
- **Score:** P=7, U=9, V=7 → **8**

### 107. hdx_back_exchange_estimator
- **Description:** Estimate per-peptide back-exchange rates from fully deuterated controls. Flag peptides with >40% back-exchange as unreliable.
- **Paper:** Same HDX-MS recommendations (DOI: 10.1038/s41592-019-0459-y)
- **pyopenms classes:** `AASequence`, `EmpiricalFormula`
- **CLI:** `--peptides peptide_list.tsv --fully-deuterated fd.tsv --max-backexchange 40 --output backexchange.tsv`
- **Domain:** proteomics (structural)
- **Score:** P=6, U=9, V=6 → **7**

### 108. xl_distance_validator
- **Description:** Given crosslinks (residue pairs) and a PDB structure, compute Cα-Cα distances. Classify as satisfied (<30Å for DSS/BS3) or violated.
- **Paper:** Distance restraints from XL-MS (DOI: 10.1002/pro.2458); XL-MS in the human cell (DOI: 10.1073/pnas.2219418120)
- **pyopenms classes:** `FASTAFile`, `AASequence`
- **CLI:** `--crosslinks crosslinks.tsv --pdb structure.pdb --crosslinker DSS --max-distance 30 --output distances.tsv`
- **Domain:** proteomics (structural)
- **Score:** P=6, U=9, V=7 → **7**

### 109. xl_link_classifier
- **Description:** Classify crosslinks as intra-protein, inter-protein, or dead-end/monolink. Compute statistics for interaction network analysis.
- **Paper:** XL-MS emerging technology (DOI: 10.1021/acs.analchem.7b04431)
- **pyopenms classes:** `FASTAFile`, `AASequence`
- **CLI:** `--crosslinks crosslinks.tsv --fasta proteome.fasta --output classified.tsv`
- **Domain:** proteomics (structural)
- **Score:** P=6, U=9, V=6 → **7**

### 110. metapeptide_lca_assigner
- **Description:** For each peptide, find all matching proteins in a metaproteome DB, retrieve taxonomy lineages, compute lowest common ancestor (LCA). Build taxonomic profile.
- **Paper:** Unipept 2024 (DOI: 10.1101/2024.09.26.615136); METATRYP v2.0 (PMID: 32897080)
- **pyopenms classes:** `FASTAFile`, `ProteaseDigestion`
- **CLI:** `--peptides identified.tsv --fasta metadb.fasta --taxonomy lineage.tsv --rank genus --output taxonomy.tsv`
- **Domain:** proteomics (metaproteomics)
- **Score:** P=7, U=9, V=7 → **8**

### 111. cleavage_site_profiler
- **Description:** From neo-N-terminal peptides (TAILS/degradomics), extract P4-P4' cleavage windows, compute position-specific amino acid frequencies for protease specificity profiling (iceLogo-style).
- **Paper:** TAILS N-terminomics (DOI: 10.1074/jbc.RA117.001113); MANTI (DOI: 10.1021/acs.analchem.1c00310)
- **pyopenms classes:** `FASTAFile`, `AASequence`, `ResidueDB`
- **CLI:** `--neo-ntermini peptides.tsv --fasta reference.fasta --window 4 --output cleavage_profile.tsv`
- **Domain:** proteomics (degradomics)
- **Score:** P=6, U=9, V=7 → **7**

### 112. nterm_modification_annotator
- **Description:** Classify N-terminal peptides as: protein N-terminus, signal peptide cleavage, transit peptide, neo-N-terminus, or acetylated. Cross-reference UniProt features.
- **Paper:** Data Processing in Positional Proteomics (DOI: 10.1002/pmic.70069)
- **pyopenms classes:** `FASTAFile`, `AASequence`, `ModificationsDB`
- **CLI:** `--input nterm_peptides.tsv --fasta reference.fasta --uniprot-features features.tsv --output annotated.tsv`
- **Domain:** proteomics (N-terminomics)
- **Score:** P=6, U=9, V=6 → **7**

### 113. immunopeptidome_qc
- **Description:** From immunopeptidome results: compute length distribution (8-12 for HLA-I), extract anchor residue frequencies (pos 2 and C-term), compute information content as QC metric.
- **Paper:** MhcVizPipe (GitHub: CaronLab/MhcVizPipe); MHC Motif Atlas (DOI: 10.1093/nar/gkac965)
- **pyopenms classes:** `AASequence`, `ResidueDB`
- **CLI:** `--input hla_peptides.tsv --hla-class I --output length_dist.tsv --motifs anchor_freq.tsv`
- **Domain:** proteomics (immunopeptidomics)
- **Score:** P=7, U=9, V=7 → **8**

### 114. proteoform_delta_annotator
- **Description:** Given observed intact proteoform masses, compute pairwise mass differences and annotate against known PTMs (acetyl +42.01, phospho +79.97, methyl +14.02, etc.).
- **Paper:** TopPIC (DOI: 10.1093/bioinformatics/btw398); ProteoCombiner (DOI: 10.1093/bioinformatics/btab175)
- **pyopenms classes:** `ModificationsDB`, `EmpiricalFormula`
- **CLI:** `--input proteoform_masses.tsv --tolerance 0.5 --output annotated_deltas.tsv`
- **Domain:** proteomics (top-down)
- **Score:** P=6, U=9, V=6 → **7**

### 115. topdown_coverage_calculator
- **Description:** From fragment ion assignments (b/y/c/z), compute per-residue bond cleavage coverage for intact proteins. Output coverage map.
- **Paper:** VisioProt-MS (DOI: 10.1177/1177932219868223)
- **pyopenms classes:** `TheoreticalSpectrumGenerator`, `AASequence`, `SpectrumAlignment`
- **CLI:** `--sequence PROTEIN_SEQ --fragments observed.tsv --ion-types b,y,c,z --tolerance 10ppm --output coverage.tsv`
- **Domain:** proteomics (top-down)
- **Score:** P=6, U=9, V=6 → **7**

### 116. library_coverage_estimator
- **Description:** Given a spectral library and FASTA proteome, compute % proteins with ≥N library entries, peptide-level coverage per protein, identify "dark" proteins with no library peptides.
- **Paper:** In silico spectral libraries (DOI: 10.1038/s41467-019-13866-z); MSLibrarian (DOI: 10.1021/acs.jproteome.1c00796)
- **pyopenms classes:** `FASTAFile`, `ProteaseDigestion`, `AASequence`
- **CLI:** `--library lib.tsv --fasta proteome.fasta --enzyme Trypsin --min-peptides 2 --output coverage.tsv`
- **Domain:** proteomics (DIA)
- **Score:** P=7, U=9, V=7 → **8**

### 117. silac_halflife_calculator
- **Description:** From pulsed-SILAC time course H/L ratios, fit exponential decay model: H/L(t) = 1-exp(-kt). Compute degradation rate k and half-life = ln(2)/k. Correct for cell division.
- **Paper:** JUMPt (DOI: 10.1021/acs.analchem.1c02309); Protein lifetimes in mouse brain (DOI: 10.1038/s41467-018-06519-0)
- **pyopenms classes:** `AASequence`, `EmpiricalFormula` (SILAC label masses: heavy Lys +8.014, heavy Arg +10.008)
- **CLI:** `--input hl_ratios.tsv --timepoints 0,6,12,24,48 --doubling-time 24 --output halflives.tsv`
- **Domain:** proteomics (turnover)
- **Score:** P=7, U=9, V=7 → **8**

### 118. scp_reporter_qc
- **Description:** Single-cell proteomics QC: compute sample-to-carrier ratio (SCR) per spectrum, blank channel contamination, reporter ion intensities across TMT channels.
- **Paper:** plexDIA (DOI: 10.1038/s41587-022-01389-w); QuantQC (DOI: 10.1021/jasms.3c00238)
- **pyopenms classes:** `MSExperiment`, `MSSpectrum`, `MzMLFile`
- **CLI:** `--input run.mzML --carrier-channel 131C --single-channels 127N,127C,128N --scr-threshold 0.1 --output qc.tsv`
- **Domain:** proteomics (single-cell)
- **Score:** P=7, U=9, V=7 → **8**

### 119. biomarker_panel_roc
- **Description:** From protein quantification with case/control labels, compute per-protein ROC curves (AUC), evaluate multi-marker panels using logistic regression.
- **Paper:** CombiROC (DOI: 10.1038/srep45477); Feature selection for biomarkers (DOI: 10.1016/j.mcpro.2021.100083)
- **pyopenms classes:** Minimal — reads quantification matrices
- **CLI:** `--input protein_quant.tsv --groups case,control --max-panel-size 5 --output panel_roc.tsv`
- **Domain:** proteomics (clinical)
- **Score:** P=7, U=8, V=7 → **7**

### 120. isobaric_purity_corrector
- **Description:** Correct TMT/iTRAQ reporter ion intensities for isotopic impurity using manufacturer-provided correction factor matrices.
- **Paper:** TMT quantification methods (DOI: 10.1021/acs.jproteome.9b00227)
- **pyopenms classes:** `MSSpectrum`, `MSExperiment`
- **CLI:** `--input quantified.tsv --label TMT16plex --purity-matrix purity.csv --output corrected.tsv`
- **Domain:** proteomics (isobaric labeling)
- **Score:** P=7, U=9, V=8 → **8**

### 121. metapeptide_function_aggregator
- **Description:** Aggregate GO/KEGG/COG functional annotations from peptide-to-protein mappings, weighted by spectral counts or intensity. Build functional profile.
- **Paper:** Unipept 2024; MetaProteomeAnalyzer (DOI: 10.1021/acs.jproteome.4c00142)
- **pyopenms classes:** `FASTAFile`
- **CLI:** `--peptides identified.tsv --fasta metadb.fasta --annotations go_terms.tsv --method spectral_count --output function.tsv`
- **Domain:** proteomics (metaproteomics)
- **Score:** P=6, U=9, V=6 → **7**

### 122. protein_completeness_matrix
- **Description:** For quantification results, compute data completeness (% non-missing values) per protein and per sample. Flag proteins below completeness threshold.
- **Paper:** SCoPE2 computational analysis (scp.slavovlab.net); DEP Bioconductor package analysis
- **pyopenms classes:** Minimal
- **CLI:** `--input quant_matrix.tsv --min-completeness 0.5 --output completeness.tsv`
- **Domain:** proteomics
- **Score:** P=6, U=8, V=7 → **7**

---

## Category 12: Niche Metabolomics — Paper-Level Tools (18 tools)

### 123. van_krevelen_data_generator
- **Description:** Compute H:C and O:C ratios from molecular formulas for Van Krevelen diagram analysis. Classify compounds into biochemical classes (lipids, amino acids, carbohydrates, nucleotides, condensed hydrocarbons).
- **Paper:** Brockman et al. "OpenVanKrevelen" (DOI: 10.1007/s11306-018-1343-y)
- **pyopenms classes:** `EmpiricalFormula`
- **CLI:** `--input formulas.tsv --classify --output van_krevelen.tsv`
- **Domain:** metabolomics
- **Score:** P=7, U=9, V=7 → **8**

### 124. kendrick_mass_defect_analyzer
- **Description:** Compute Kendrick mass defect for arbitrary repeating units (CH2, CF2, C2H4O). Group features into homologous series. Detect lipid chains, PFAS, polymers.
- **Paper:** Sleno "Mass defect in modern MS" (J Mass Spectrom 2012); Fouquet & Sato (2017)
- **pyopenms classes:** `EmpiricalFormula`
- **CLI:** `--input features.tsv --base CH2 --output kmd_series.tsv`
- **Domain:** metabolomics
- **Score:** P=7, U=9, V=7 → **8**

### 125. drug_metabolite_screener
- **Description:** Given a parent drug formula, generate expected metabolite masses from Phase I (oxidation +15.995, demethylation -14.016, hydrolysis +18.011) and Phase II (glucuronidation +176.032, sulfation +79.957, glutathione +305.068) reactions. Screen mzML against predictions.
- **Paper:** Zhang et al. "Mass defect filter for drug metabolites" (PMID: 19598168); BioTransformer (DOI: 10.1186/s13321-018-0324-5)
- **pyopenms classes:** `EmpiricalFormula`, `MSExperiment`, `MzMLFile`
- **CLI:** `--parent-formula C17H14ClN3O --reactions phase1,phase2 --input run.mzML --ppm 5 --output metabolites.tsv`
- **Domain:** metabolomics (drug metabolism)
- **Score:** P=7, U=9, V=8 → **8**

### 126. isf_detector
- **Description:** Detect in-source fragmentation artifacts by coelution analysis (RT correlation ≥0.9) and matching fragment ions in MS2 of putative precursors. Common neutral losses: H2O, CO2, NH3.
- **Paper:** ISFrag (DOI: 10.1021/acs.analchem.1c01644); Mahieu et al. (DOI: 10.1021/ac504118y)
- **pyopenms classes:** `MSExperiment`, `MzMLFile`, `MSSpectrum`
- **CLI:** `--input features.tsv --mzml run.mzML --rt-correlation 0.9 --output isf_annotated.tsv`
- **Domain:** metabolomics
- **Score:** P=7, U=9, V=8 → **8**

### 127. mid_natural_abundance_corrector
- **Description:** Correct mass isotopomer distributions (MID) for natural 13C/15N/2H abundance. Build correction matrix from molecular formula, solve via NNLS. Essential for 13C flux analysis.
- **Paper:** Corna (DOI: 10.1101/2020.09.19.304741); Fernandez et al. (J Mass Spectrom 1996)
- **pyopenms classes:** `EmpiricalFormula`, `CoarseIsotopePatternGenerator`
- **CLI:** `--input isotopologues.tsv --formula C6H12O6 --tracer 13C --method nnls --output corrected_mid.tsv`
- **Domain:** metabolomics (fluxomics)
- **Score:** P=7, U=9, V=8 → **8**

### 128. isotope_label_detector
- **Description:** Compare unlabeled vs. uniformly 13C/15N-labeled samples. Detect paired features by RT and expected mass shift. Report carbon/nitrogen count per metabolite.
- **Paper:** X13CMS (DOI: 10.1021/ac403384n); MetExtract II (DOI: 10.1021/acs.analchem.7b02518)
- **pyopenms classes:** `MzMLFile`, `EmpiricalFormula`, `CoarseIsotopePatternGenerator`
- **CLI:** `--unlabeled unlabeled.mzML --labeled labeled.mzML --tracer 13C --ppm 5 --output labeled_pairs.tsv`
- **Domain:** metabolomics (stable isotope tracing)
- **Score:** P=7, U=9, V=7 → **8**

### 129. formula_validator_golden_rules
- **Description:** Apply Kind & Fiehn's Seven Golden Rules to filter molecular formula candidates: element limits, LEWIS/SENIOR valence, H/C ratio (0.2-3.1), N/C O/C S/C P/C ratios, RDBE check.
- **Paper:** Kind & Fiehn "Seven Golden Rules" (DOI: 10.1186/1471-2105-8-105) — cited >2500 times
- **pyopenms classes:** `EmpiricalFormula`
- **CLI:** `--input candidate_formulas.tsv --rules all --output validated.tsv`
- **Domain:** metabolomics
- **Score:** P=8, U=9, V=8 → **8**

### 130. rdbe_calculator
- **Description:** Calculate Ring/Double Bond Equivalence (RDBE = (2C+2-H+N+P)/2) for molecular formulas. Flag impossible formulas (negative RDBE, non-integer for even-electron ions).
- **Paper:** Fiehn Lab documentation; universally used in structure elucidation
- **pyopenms classes:** `EmpiricalFormula`
- **CLI:** `--input formulas.tsv --output rdbe_values.tsv`
- **Domain:** metabolomics
- **Score:** P=6, U=8, V=7 → **7**

### 131. spectral_entropy_scorer
- **Description:** Weight MS/MS peaks by Shannon entropy, compute entropy-based similarity instead of cosine. Outperforms 42 alternative metrics per the Li & Fiehn paper.
- **Paper:** Li & Fiehn (DOI: 10.1038/s41592-021-01331-z) — Nature Methods 2021
- **pyopenms classes:** `MSSpectrum`, `MzMLFile`
- **CLI:** `--query query.mgf --library ref.mgf --tolerance 0.02 --output entropy_scores.tsv`
- **Domain:** metabolomics / proteomics
- **Score:** P=8, U=9, V=8 → **8**

### 132. suspect_screener
- **Description:** Match detected exact masses against curated suspect screening lists (CompTox, NORMAN, PubChemLite). Rank by mass accuracy, isotope fit, data source count.
- **Paper:** Scannotation (DOI: 10.1021/acs.est.3c04764); EPA CompTox Dashboard
- **pyopenms classes:** `EmpiricalFormula`
- **CLI:** `--input features.tsv --suspects comptox_list.csv --ppm 5 --output matches.tsv`
- **Domain:** metabolomics (exposomics/environmental)
- **Score:** P=7, U=9, V=7 → **8**

### 133. lipid_ecn_rt_predictor
- **Description:** Predict lipid RT from Equivalent Carbon Number (ECN = carbons - 2×double bonds). Lipids with same headgroup show linear RT vs ECN. Transfer predictions across lipid classes.
- **Paper:** Koelmel et al. (DOI: 10.1021/acs.analchem.1c03770)
- **pyopenms classes:** `EmpiricalFormula`
- **CLI:** `--input lipid_annotations.tsv --calibration standards.tsv --output rt_predictions.tsv`
- **Domain:** metabolomics (lipidomics)
- **Score:** P=6, U=9, V=6 → **7**

### 134. lipid_species_resolver
- **Description:** From sum composition lipid annotations (e.g., "PC 36:2"), enumerate all acyl chain combinations (16:0/20:2, 18:1/18:1), compute masses, match against observed MS2 fragments.
- **Paper:** Goslin nomenclature (DOI: 10.1021/acs.analchem.0c01690); Fatty acyl C=C (DOI: 10.1038/s41467-025-61911-x)
- **pyopenms classes:** `EmpiricalFormula`, `TheoreticalSpectrumGenerator` concepts
- **CLI:** `--input lipid_features.tsv --lipid-class PC --fragments fragments.tsv --ppm 5 --output resolved.tsv`
- **Domain:** metabolomics (lipidomics)
- **Score:** P=6, U=9, V=6 → **7**

### 135. kovats_ri_calculator
- **Description:** Calculate Kovats Retention Index from alkane standard RTs. Convert instrument-specific RT to universal RI for NIST/Fiehnlib library matching.
- **Paper:** Kind et al. (DOI: 10.1093/bioinformatics/btp056)
- **pyopenms classes:** `MzMLFile`
- **CLI:** `--input features.tsv --alkane-standards standards.tsv --mode temperature_programmed --output ri_converted.tsv`
- **Domain:** metabolomics (GC-MS)
- **Score:** P=6, U=9, V=6 → **7**

### 136. gnps_fbmn_exporter
- **Description:** Export consensus MS2 spectra + quantification table in GNPS Feature-Based Molecular Networking format. Handle SCANS, PEPMASS, RTINSECONDS, MSLEVEL fields correctly.
- **Paper:** Nothias et al. FBMN (DOI: 10.1038/s41592-020-0933-6)
- **pyopenms classes:** `MzMLFile`, `FeatureXMLFile`, `MSSpectrum`
- **CLI:** `--mzml-dir runs/ --features features.tsv --output-mgf gnps.mgf --output-quant gnps_quant.csv`
- **Domain:** metabolomics (molecular networking)
- **Score:** P=8, U=8, V=8 → **8**

### 137. isotope_pattern_fit_scorer
- **Description:** Compare experimental isotope envelope (M, M+1, M+2, M+3) against theoretical distribution for candidate formulas. Score by chi-square, cosine, sigma. Detect Cl/Br/S signatures from M+2 intensity.
- **Paper:** Pluskal et al. (DOI: 10.1021/ac3000418); Kind & Fiehn Seven Golden Rules (Rule 3)
- **pyopenms classes:** `EmpiricalFormula`, `CoarseIsotopePatternGenerator`, `FineIsotopePatternGenerator`
- **CLI:** `--observed "180.063:100,181.067:6.5,182.070:0.5" --candidates formulas.tsv --charge 1 --output scored.tsv`
- **Domain:** metabolomics
- **Score:** P=7, U=8, V=8 → **8**

### 138. mass_difference_network_builder
- **Description:** Build mass difference network connecting features whose Δm/z matches known biotransformations. Detect oxidation (+15.995), methylation (+14.016), glucuronidation (+176.032), etc.
- **Paper:** GNPS molecular networking principle; Watrous et al. (DOI: 10.1073/pnas.1203689109)
- **pyopenms classes:** `EmpiricalFormula`
- **CLI:** `--input features.tsv --reactions biotransformations.tsv --tolerance 0.005 --output network.tsv`
- **Domain:** metabolomics
- **Score:** P=7, U=9, V=7 → **8**

### 139. metabolite_class_predictor
- **Description:** Predict compound class (lipid, amino acid, nucleotide, carbohydrate, organic acid) from mass defect ranges, H:C/O:C ratios, and RDBE values — without structural information.
- **Paper:** ClassyFire (DOI: 10.1186/s13321-016-0174-y) uses structure; this uses mass-only heuristics
- **pyopenms classes:** `EmpiricalFormula`
- **CLI:** `--input formulas.tsv --output class_predictions.tsv`
- **Domain:** metabolomics
- **Score:** P=6, U=9, V=6 → **7**

### 140. blank_subtraction_tool
- **Description:** Compare sample vs. blank features by m/z+RT matching. Remove features present in blank or below fold-change threshold. Standard metabolomics preprocessing.
- **Paper:** UmetaFlow (DOI: 10.1186/s13321-023-00724-w); XCMS-based workflows
- **pyopenms classes:** `MSExperiment`, `MzMLFile`
- **CLI:** `--sample sample_features.tsv --blank blank_features.tsv --fold-change 3 --mz-ppm 5 --rt-tolerance 10 --output cleaned.tsv`
- **Domain:** metabolomics
- **Score:** P=7, U=9, V=8 → **8**

---

## Category 13: Enhanced Features for Existing Tools

*Paper-level feature additions that enrich tools in Categories 1-10:*

### Tool #2 spectrum_similarity_scorer — Enhanced
- Add **spectral entropy scoring** (Li & Fiehn 2021, Nature Methods) alongside cosine
- Add **modified cosine** for shifted precursors (molecular networking)
- Add **Stein-Scott composite score** for library matching
- Add **neutral loss-aware matching** for metabolomics

### Tool #9 xic_extractor — Enhanced
- Add **peak area integration** (trapezoidal rule) — reported in every quantification paper
- Add **FWHM calculation** per XIC peak — standard chromatographic QC metric
- Add **signal-to-noise estimation** — essential for LOD/LOQ determination
- Add **batch mode** for extracting XICs across multiple runs simultaneously

### Tool #16 peptide_property_calculator — Enhanced
- Add **Kyte-Doolittle sliding window** hydrophobicity profiles — amphipathicity analysis
- Add **Boman interaction index** — antimicrobial peptide studies (modlAMP feature)
- Add **Wimley-White interfacial hydrophobicity** — membrane interaction prediction
- Add **net charge vs pH curve** generation — electrophoresis prediction

### Tool #22 missed_cleavage_analyzer — Enhanced
- Add **P4-P4' position-specific frequency analysis** — iceLogo-style cleavage specificity
- Add **enzyme specificity profiling** — beyond missed cleavages, analyze actual vs. expected cleavage sites
- Add **comparison mode** — compare missed cleavage rates across multiple runs

### Tool #50 lc_ms_qc_reporter — Enhanced
- Add **longitudinal QC trending** — track metrics across a batch of runs over time
- Add **system suitability criteria** — pass/fail thresholds for each metric
- Add **peak capacity estimation** from MS1 data — method optimization metric
- Add **TopN sampling depth** — fraction of MS1 features with MS2 attempts

### Tool #58 adduct_calculator — Enhanced
- Add **in-source fragmentation prediction** — common ISF losses alongside adducts
- Add **multimer prediction** — [2M+H]+, [2M+Na]+, [3M+H]+
- Add **cluster ion prediction** — [M+Na+K-H]+, [M+2Na-H]+

---

## Score Summary (All 145 Tools)

| Score | Count | Examples |
|-------|-------|---------|
| **9**  | 5 | spectrum_similarity_scorer, xic_extractor, peptide_property_calculator, fasta_subset_extractor, adduct_calculator |
| **8**  | 54 | + niche tools: phosphosite_class_filter, hdx_deuterium_uptake, golden_rules_validator, spectral_entropy, drug_metabolite_screener, FBMN exporter, MID corrector, etc. |
| **7**  | 58 | + niche tools: xl_distance_validator, cleavage_site_profiler, lipid_species_resolver, kovats_ri, etc. |
| **6**  | 18 | score-6 utilities and emerging tools |

## Implementation Priority Tiers

**Tier 1 — Implement First (Score 9, highest impact — 5 tools):**
1. spectrum_similarity_scorer (#2)
2. xic_extractor (#9)
3. peptide_property_calculator (#16)
4. fasta_subset_extractor (#32)
5. adduct_calculator (#58)

**Tier 2 — High Value (Score 8 — 54 tools):**
Core utilities: #1, #3, #5, #7, #8, #11, #14, #15, #17, #18, #19, #22, #24, #25, #26, #33, #34, #40, #42, #43, #48, #49, #50, #51, #52, #54, #56, #59, #62, #64, #73, #82, #83, #84, #85
Niche paper tools: #101-103, #104, #106, #110, #113, #116-118, #120, #123-129, #131, #132, #136-138, #140

**Tier 3 — Solid Utilities (Score 7 — 58 tools):**
Core: #4, #6, #10, #12, #13, #20, #21, #23, #27-29, #35-38, #41, #44, #47, #53, #55, #57, #60, #61, #63, #65-69, #74, #75, #78, #80, #81, #86-91, #93-97, #100
Niche: #105, #107-109, #111, #112, #114, #115, #119, #121, #122, #130, #133-135, #139

**Tier 4 — Nice to Have (Score 6 — 18 tools):**
#30, #31, #39, #45, #46, #71, #72, #76, #77, #79, #92, #98, #99 + remaining niche
