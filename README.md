# Agentomics

A growing collection of **123 standalone CLI tools** built with [pyopenms](https://pyopenms.readthedocs.io/) for proteomics and metabolomics workflows. Every tool in this repository fills a gap not covered by existing OpenMS TOPP tools — small, focused utilities that researchers need daily but typically write as throwaway scripts.

## Why This Exists

Mass spectrometry researchers constantly need small utilities: extract an XIC from an mzML file, compute adduct m/z values for a metabolite, check peptide uniqueness in a FASTA database, validate crosslink distances against a PDB structure. These tasks are too simple for a full pipeline but too tedious to re-implement from scratch every time.

Agentomics collects these utilities into a single, organized repository where each tool is:

- **Self-contained** — no cross-tool dependencies, install and run independently
- **CLI-first** — every tool has an `click` interface, usable from the command line or imported as a Python library
- **Tested** — every tool ships with unit tests using synthetic pyopenms data
- **pyopenms-native** — built on the official Python bindings for OpenMS, not reimplementing what already exists

## AI-Generated Disclaimer

> **All code in this repository is written entirely by AI agents** (Claude Code, GitHub Copilot, Cursor, Gemini, etc.). This is an agentic-only development project — tool ideas were researched from GitHub repositories, community forums (BioStars, Reddit), published papers, and pyopenms documentation, then implemented by AI. Human review is applied for quality control and direction, but the code itself is machine-generated. Use at your own discretion and always validate results against established tools for critical analyses.

## Contributing (Agentic Workflow)

This repo is designed for AI agent contributions. The full contributor guide is in [AGENTS.md](AGENTS.md), but the key idea is:

1. **Pick a gap** — find a utility task that researchers need but no TOPP tool covers
2. **Follow the structure** — every tool lives in its own directory with a standard layout (see below)
3. **Validate in isolation** — each tool must pass `ruff check` and `pytest` in a fresh venv with only `pyopenms` installed
4. **Do not duplicate TOPP tools** — if `FileConverter`, `PeakPickerHiRes`, `FalseDiscoveryRate`, or any other TOPP command already does it, don't rebuild it here

Two Claude Code skills are available for contributors:
- **`contribute-script`** — guided workflow for adding a new tool
- **`validate-script`** — validate any tool in an isolated venv (ruff + pytest)

## Tool Structure

Every tool follows the same directory layout:

```
tools/<domain>/<topic>/<tool_name>/
├── <tool_name>.py        # The tool (importable functions + click CLI)
├── requirements.txt      # pyopenms + tool-specific deps (no version pins)
├── README.md             # Brief description + CLI usage examples
└── tests/
    ├── conftest.py       # requires_pyopenms marker + sys.path setup
    └── test_<tool_name>.py
```

**Every `.py` file contains:**

1. A module docstring describing the tool, its features, and usage
2. A pyopenms import guard with a user-friendly error message
3. Importable functions with type hints and numpy-style docstrings — so the tool works both as a library and as a CLI
4. A `main()` function wiring up `click` for command-line usage
5. An `if __name__ == "__main__": main()` guard

**Domains:** `proteomics/`, `metabolomics/`

**Proteomics topics:** `spectrum_analysis/`, `peptide_analysis/`, `protein_analysis/`, `fasta_utils/`, `file_conversion/`, `quality_control/`, `targeted_proteomics/`, `identification/`, `ptm_analysis/`, `structural_proteomics/`, `specialized/`, `rna/`

**Metabolomics topics:** `formula_tools/`, `feature_processing/`, `spectral_analysis/`, `compound_annotation/`, `drug_metabolism/`, `isotope_labeling/`, `lipidomics/`, `export/`

## Requirements

```bash
pip install pyopenms
```

Some tools require additional dependencies (`numpy`, `scipy`). Check each tool's `requirements.txt`.

## Running a Tool

```bash
# Install dependencies
pip install -r tools/proteomics/spectrum_analysis/theoretical_spectrum_generator/requirements.txt

# Run via CLI
python tools/proteomics/spectrum_analysis/theoretical_spectrum_generator/theoretical_spectrum_generator.py --help

# Run tests
PYTHONPATH=tools/proteomics/spectrum_analysis/theoretical_spectrum_generator \
  python -m pytest tools/proteomics/spectrum_analysis/theoretical_spectrum_generator/tests/ -v
```

## Validation

Each tool is validated in an isolated venv:

```bash
TOOL_DIR=tools/<domain>/<topic>/<tool_name>
VENV_DIR=$(mktemp -d)
python -m venv "$VENV_DIR"
"$VENV_DIR/bin/python" -m pip install -r "$TOOL_DIR/requirements.txt"
"$VENV_DIR/bin/python" -m pip install pytest ruff
"$VENV_DIR/bin/python" -m ruff check "$TOOL_DIR/"
PYTHONPATH="$TOOL_DIR" "$VENV_DIR/bin/python" -m pytest "$TOOL_DIR/tests/" -v
rm -rf "$VENV_DIR"
```

Both `ruff` and `pytest` must pass with zero errors.

---

## Tool Catalog

### Proteomics (89 tools)

#### Spectrum Analysis (7 tools)

| Tool | Description |
|------|-------------|
| [`theoretical_spectrum_generator`](tools/proteomics/spectrum_analysis/theoretical_spectrum_generator/) | Generate theoretical b/y/a/c/x/z fragment ion spectra for peptide sequences |
| [`spectrum_similarity_scorer`](tools/proteomics/spectrum_analysis/spectrum_similarity_scorer/) | Compute cosine similarity between MS2 spectra from MGF files |
| [`spectrum_annotator`](tools/proteomics/spectrum_analysis/spectrum_annotator/) | Annotate observed MS2 peaks with theoretical fragment ion matches |
| [`spectrum_scoring_hyperscore`](tools/proteomics/spectrum_analysis/spectrum_scoring_hyperscore/) | Score experimental spectra against theoretical using HyperScore |
| [`spectrum_entropy_calculator`](tools/proteomics/spectrum_analysis/spectrum_entropy_calculator/) | Calculate normalized Shannon entropy for MS2 spectra |
| [`spectral_library_builder`](tools/proteomics/spectrum_analysis/spectral_library_builder/) | Build consensus spectral libraries from mzML + peptide identifications |
| [`spectral_library_format_converter`](tools/proteomics/spectrum_analysis/spectral_library_format_converter/) | Convert between spectral library formats (MSP, TraML) |

#### Peptide Analysis (12 tools)

| Tool | Description |
|------|-------------|
| [`peptide_property_calculator`](tools/proteomics/peptide_analysis/peptide_property_calculator/) | Calculate pI, hydrophobicity, charge at pH, amino acid composition |
| [`peptide_mass_calculator`](tools/proteomics/peptide_analysis/peptide_mass_calculator/) | Monoisotopic/average masses and b/y fragment ions |
| [`peptide_uniqueness_checker`](tools/proteomics/peptide_analysis/peptide_uniqueness_checker/) | Check if peptides are proteotypic within a FASTA database |
| [`modification_mass_calculator`](tools/proteomics/peptide_analysis/modification_mass_calculator/) | Query Unimod by name or mass shift, compute modified peptide masses |
| [`modified_peptide_generator`](tools/proteomics/peptide_analysis/modified_peptide_generator/) | Enumerate all modified peptide variants for given variable/fixed mods |
| [`peptide_modification_analyzer`](tools/proteomics/peptide_analysis/peptide_modification_analyzer/) | Residue-by-residue mass breakdown of modified peptides |
| [`peptide_detectability_predictor`](tools/proteomics/peptide_analysis/peptide_detectability_predictor/) | Predict peptide detectability from physicochemical heuristics |
| [`isoelectric_point_calculator`](tools/proteomics/peptide_analysis/isoelectric_point_calculator/) | Calculate pI using Henderson-Hasselbalch with configurable pK sets |
| [`charge_state_predictor`](tools/proteomics/peptide_analysis/charge_state_predictor/) | Predict charge state distribution based on basic residues |
| [`amino_acid_composition_analyzer`](tools/proteomics/peptide_analysis/amino_acid_composition_analyzer/) | Amino acid frequency and composition statistics |
| [`rt_prediction_additive`](tools/proteomics/peptide_analysis/rt_prediction_additive/) | Predict peptide RT using additive hydrophobicity models |
| [`peptide_mass_fingerprint`](tools/proteomics/peptide_analysis/peptide_mass_fingerprint/) | Generate/match peptide mass fingerprints for MALDI-TOF identification |

#### Protein Analysis (5 tools)

| Tool | Description |
|------|-------------|
| [`protein_digest`](tools/proteomics/protein_analysis/protein_digest/) | In-silico enzymatic protein digestion |
| [`protein_coverage_calculator`](tools/proteomics/protein_analysis/protein_coverage_calculator/) | Map peptides to proteins and calculate sequence coverage |
| [`protein_group_reporter`](tools/proteomics/protein_analysis/protein_group_reporter/) | Generate clean protein-level reports with group membership |
| [`spectral_counting_quantifier`](tools/proteomics/protein_analysis/spectral_counting_quantifier/) | Calculate protein abundances using emPAI or NSAF methods |
| [`peptide_to_protein_mapper`](tools/proteomics/protein_analysis/peptide_to_protein_mapper/) | Map peptide sequences to parent proteins in a FASTA database |

#### FASTA Utilities (8 tools)

| Tool | Description |
|------|-------------|
| [`fasta_subset_extractor`](tools/proteomics/fasta_utils/fasta_subset_extractor/) | Extract proteins by accession list, keyword, or length range |
| [`fasta_statistics_reporter`](tools/proteomics/fasta_utils/fasta_statistics_reporter/) | Report protein count, lengths, amino acid frequency, tryptic peptide counts |
| [`contaminant_database_merger`](tools/proteomics/fasta_utils/contaminant_database_merger/) | Append cRAP contaminant sequences with configurable prefix |
| [`fasta_cleaner`](tools/proteomics/fasta_utils/fasta_cleaner/) | Remove duplicates, fix headers, filter by length |
| [`fasta_merger`](tools/proteomics/fasta_utils/fasta_merger/) | Merge multiple FASTA files with duplicate removal |
| [`fasta_decoy_validator`](tools/proteomics/fasta_utils/fasta_decoy_validator/) | Check if a FASTA already contains decoys, validate prefix consistency |
| [`fasta_in_silico_digest_stats`](tools/proteomics/fasta_utils/fasta_in_silico_digest_stats/) | Digest a FASTA and report peptide-level statistics |
| [`fasta_taxonomy_splitter`](tools/proteomics/fasta_utils/fasta_taxonomy_splitter/) | Split multi-organism FASTA by taxonomy from headers |

#### File Conversion (8 tools)

| Tool | Description |
|------|-------------|
| [`mzml_to_mgf_converter`](tools/proteomics/file_conversion/mzml_to_mgf_converter/) | Convert MS2 spectra from mzML to MGF format |
| [`mgf_to_mzml_converter`](tools/proteomics/file_conversion/mgf_to_mzml_converter/) | Convert MGF files to mzML format |
| [`consensus_map_to_matrix`](tools/proteomics/file_conversion/consensus_map_to_matrix/) | Convert consensusXML to flat quantification matrix |
| [`idxml_to_tsv_exporter`](tools/proteomics/file_conversion/idxml_to_tsv_exporter/) | Export idXML identification results to flat TSV |
| [`ms_data_to_csv_exporter`](tools/proteomics/file_conversion/ms_data_to_csv_exporter/) | Export mzML/featureXML data to CSV with column selection |
| [`mztab_summarizer`](tools/proteomics/file_conversion/mztab_summarizer/) | Parse mzTab files and extract summary statistics |
| [`featurexml_merger`](tools/proteomics/file_conversion/featurexml_merger/) | Merge multiple featureXML files |
| [`ms_data_ml_exporter`](tools/proteomics/file_conversion/ms_data_ml_exporter/) | Export MS features as ML-ready matrices |

#### Quality Control (15 tools)

| Tool | Description |
|------|-------------|
| [`lc_ms_qc_reporter`](tools/proteomics/quality_control/lc_ms_qc_reporter/) | Comprehensive QC report from mzML (TIC, MS1/MS2 counts, charge distribution) |
| [`mzqc_generator`](tools/proteomics/quality_control/mzqc_generator/) | Generate mzQC-format (HUPO-PSI standard) quality control files |
| [`identification_qc_reporter`](tools/proteomics/quality_control/identification_qc_reporter/) | Report identification-level QC metrics from search results |
| [`run_comparison_reporter`](tools/proteomics/quality_control/run_comparison_reporter/) | Compare mzML files side-by-side (TIC correlation, shared precursors) |
| [`mass_error_distribution_analyzer`](tools/proteomics/quality_control/mass_error_distribution_analyzer/) | Compute precursor and fragment mass error distributions |
| [`acquisition_rate_analyzer`](tools/proteomics/quality_control/acquisition_rate_analyzer/) | Analyze MS1/MS2 acquisition rates, cycle time, duty cycle |
| [`precursor_isolation_purity`](tools/proteomics/quality_control/precursor_isolation_purity/) | Estimate precursor isolation purity and co-isolation interference |
| [`injection_time_analyzer`](tools/proteomics/quality_control/injection_time_analyzer/) | Extract and analyze injection time values from mzML |
| [`collision_energy_analyzer`](tools/proteomics/quality_control/collision_energy_analyzer/) | Extract and analyze collision energy values across MS2 spectra |
| [`precursor_charge_distribution`](tools/proteomics/quality_control/precursor_charge_distribution/) | Analyze charge state distribution across MS2 spectra |
| [`precursor_recurrence_analyzer`](tools/proteomics/quality_control/precursor_recurrence_analyzer/) | Analyze precursor resampling frequency in DDA runs |
| [`missed_cleavage_analyzer`](tools/proteomics/quality_control/missed_cleavage_analyzer/) | Analyze missed cleavage distribution as a digestion QC metric |
| [`sample_complexity_estimator`](tools/proteomics/quality_control/sample_complexity_estimator/) | Estimate sample complexity from MS1 peak density |
| [`spectrum_file_info`](tools/proteomics/quality_control/spectrum_file_info/) | Summary statistics for mzML files |
| [`ms1_feature_intensity_tracker`](tools/proteomics/quality_control/ms1_feature_intensity_tracker/) | Track feature intensities across a batch of mzML runs |

#### Targeted Proteomics (7 tools)

| Tool | Description |
|------|-------------|
| [`xic_extractor`](tools/proteomics/targeted_proteomics/xic_extractor/) | Extract ion chromatograms for target m/z values from mzML |
| [`tic_bpc_calculator`](tools/proteomics/targeted_proteomics/tic_bpc_calculator/) | Compute TIC and base peak chromatograms from mzML |
| [`transition_list_generator`](tools/proteomics/targeted_proteomics/transition_list_generator/) | Generate SRM/MRM/PRM transition lists from peptide sequences |
| [`irt_calculator`](tools/proteomics/targeted_proteomics/irt_calculator/) | Convert observed RT to indexed retention time (iRT) values |
| [`inclusion_list_generator`](tools/proteomics/targeted_proteomics/inclusion_list_generator/) | Generate instrument inclusion lists from identification results |
| [`dia_window_analyzer`](tools/proteomics/targeted_proteomics/dia_window_analyzer/) | Report DIA isolation window scheme from mzML metadata |
| [`library_coverage_estimator`](tools/proteomics/targeted_proteomics/library_coverage_estimator/) | Estimate proteome coverage of a spectral library |

#### Identification (7 tools)

| Tool | Description |
|------|-------------|
| [`feature_detection_proteomics`](tools/proteomics/identification/feature_detection_proteomics/) | Peptide feature detection from LC-MS/MS data |
| [`psm_feature_extractor`](tools/proteomics/identification/psm_feature_extractor/) | Extract rescoring features from PSMs (mass error, coverage, intensity) |
| [`peptide_spectral_match_validator`](tools/proteomics/identification/peptide_spectral_match_validator/) | Validate individual PSMs by recomputing fragment ion coverage |
| [`semi_tryptic_peptide_finder`](tools/proteomics/identification/semi_tryptic_peptide_finder/) | Classify peptides as fully/semi/non-tryptic |
| [`sequence_tag_generator`](tools/proteomics/identification/sequence_tag_generator/) | Generate de novo sequence tags from MS2 fragment ion ladders |
| [`mzml_spectrum_subsetter`](tools/proteomics/identification/mzml_spectrum_subsetter/) | Extract specific spectra from mzML by scan number list |
| [`mzml_metadata_extractor`](tools/proteomics/identification/mzml_metadata_extractor/) | Extract instrument metadata from mzML files |

#### PTM Analysis (5 tools)

| Tool | Description |
|------|-------------|
| [`ptm_site_localization_scorer`](tools/proteomics/ptm_analysis/ptm_site_localization_scorer/) | Score PTM site localization confidence using fragment ion coverage |
| [`phosphosite_class_filter`](tools/proteomics/ptm_analysis/phosphosite_class_filter/) | Classify phosphosites into Class I/II/III by localization probability |
| [`phospho_motif_analyzer`](tools/proteomics/ptm_analysis/phospho_motif_analyzer/) | Extract sequence windows around phosphosites and analyze kinase motifs |
| [`phospho_enrichment_qc`](tools/proteomics/ptm_analysis/phospho_enrichment_qc/) | Compute phospho-enrichment efficiency and pSer/pThr/pTyr ratios |
| [`glycopeptide_mass_calculator`](tools/proteomics/ptm_analysis/glycopeptide_mass_calculator/) | Calculate glycopeptide masses with glycan compositions |

#### Structural Proteomics (5 tools)

| Tool | Description |
|------|-------------|
| [`hdx_deuterium_uptake`](tools/proteomics/structural_proteomics/hdx_deuterium_uptake/) | Calculate deuterium uptake from HDX-MS time course data |
| [`hdx_back_exchange_estimator`](tools/proteomics/structural_proteomics/hdx_back_exchange_estimator/) | Estimate per-peptide back-exchange rates from fully deuterated controls |
| [`crosslink_mass_calculator`](tools/proteomics/structural_proteomics/crosslink_mass_calculator/) | Calculate masses for crosslinked peptide pairs (DSS, BS3, DSSO) |
| [`xl_distance_validator`](tools/proteomics/structural_proteomics/xl_distance_validator/) | Validate crosslink distances against PDB structures |
| [`xl_link_classifier`](tools/proteomics/structural_proteomics/xl_link_classifier/) | Classify crosslinks as intra-protein, inter-protein, or monolink |

#### Specialized (7 tools)

| Tool | Description |
|------|-------------|
| [`immunopeptide_filter`](tools/proteomics/specialized/immunopeptide_filter/) | Filter peptides for MHC-I/II by length range and motif |
| [`immunopeptidome_qc`](tools/proteomics/specialized/immunopeptidome_qc/) | QC for immunopeptidomics (length distribution, anchor residues) |
| [`metapeptide_lca_assigner`](tools/proteomics/specialized/metapeptide_lca_assigner/) | Assign lowest common ancestor taxonomy from peptide-protein mappings |
| [`cleavage_site_profiler`](tools/proteomics/specialized/cleavage_site_profiler/) | Profile protease cleavage site specificity from N-terminomics data |
| [`nterm_modification_annotator`](tools/proteomics/specialized/nterm_modification_annotator/) | Classify N-terminal peptides (protein N-term, signal peptide, neo-N-term) |
| [`proteoform_delta_annotator`](tools/proteomics/specialized/proteoform_delta_annotator/) | Annotate mass differences between proteoforms with known PTMs |
| [`topdown_coverage_calculator`](tools/proteomics/specialized/topdown_coverage_calculator/) | Compute per-residue bond cleavage coverage for intact proteins |

#### RNA (3 tools)

| Tool | Description |
|------|-------------|
| [`rna_mass_calculator`](tools/proteomics/rna/rna_mass_calculator/) | Calculate mass, formula, and isotopes for RNA/oligonucleotide sequences |
| [`rna_digest`](tools/proteomics/rna/rna_digest/) | In silico RNA digestion with RNases (T1, U2, etc.) |
| [`rna_fragment_spectrum_generator`](tools/proteomics/rna/rna_fragment_spectrum_generator/) | Generate theoretical RNA fragment spectra (c/y/w/a-B ions) |

---

### Metabolomics (34 tools)

#### Formula Tools (8 tools)

| Tool | Description |
|------|-------------|
| [`adduct_calculator`](tools/metabolomics/formula_tools/adduct_calculator/) | Compute m/z for all common ESI adducts given a formula or mass |
| [`molecular_formula_finder`](tools/metabolomics/formula_tools/molecular_formula_finder/) | Enumerate valid molecular formulas for an accurate mass with element constraints |
| [`mass_decomposition_tool`](tools/metabolomics/formula_tools/mass_decomposition_tool/) | Find molecular formula compositions for a given mass within tolerance |
| [`formula_mass_calculator`](tools/metabolomics/formula_tools/formula_mass_calculator/) | Calculate exact masses for molecular formulas with adduct support |
| [`formula_validator_golden_rules`](tools/metabolomics/formula_tools/formula_validator_golden_rules/) | Apply Kind & Fiehn's Seven Golden Rules to filter formula candidates |
| [`rdbe_calculator`](tools/metabolomics/formula_tools/rdbe_calculator/) | Calculate Ring/Double Bond Equivalence for molecular formulas |
| [`metabolite_formula_annotator`](tools/metabolomics/formula_tools/metabolite_formula_annotator/) | Annotate features with candidate formulas using mass + isotope fit |
| [`mass_accuracy_calculator`](tools/metabolomics/formula_tools/mass_accuracy_calculator/) | Compute m/z mass accuracy (ppm error) for sequences or formulas |

#### Feature Processing (7 tools)

| Tool | Description |
|------|-------------|
| [`blank_subtraction_tool`](tools/metabolomics/feature_processing/blank_subtraction_tool/) | Subtract blank/control features from sample features by m/z + RT matching |
| [`duplicate_feature_detector`](tools/metabolomics/feature_processing/duplicate_feature_detector/) | Detect and flag duplicate features by m/z and RT proximity |
| [`adduct_group_analyzer`](tools/metabolomics/feature_processing/adduct_group_analyzer/) | Group features by adduct relationships into ion identity groups |
| [`isf_detector`](tools/metabolomics/feature_processing/isf_detector/) | Detect in-source fragmentation artifacts by coelution and neutral loss |
| [`targeted_feature_extractor`](tools/metabolomics/feature_processing/targeted_feature_extractor/) | Extract features for known compounds from MS1 data |
| [`mass_defect_filter`](tools/metabolomics/feature_processing/mass_defect_filter/) | Filter features by mass defect and Kendrick mass defect |
| [`metabolite_feature_detection`](tools/metabolomics/feature_processing/metabolite_feature_detection/) | Metabolite feature detection from LC-MS data |

#### Spectral Analysis (6 tools)

| Tool | Description |
|------|-------------|
| [`spectral_entropy_scorer`](tools/metabolomics/spectral_analysis/spectral_entropy_scorer/) | Compute spectral entropy similarity (Li & Fiehn 2021) |
| [`neutral_loss_scanner`](tools/metabolomics/spectral_analysis/neutral_loss_scanner/) | Scan MS2 spectra for characteristic neutral losses |
| [`isotope_pattern_scorer`](tools/metabolomics/spectral_analysis/isotope_pattern_scorer/) | Score observed vs. theoretical isotope patterns |
| [`isotope_pattern_matcher`](tools/metabolomics/spectral_analysis/isotope_pattern_matcher/) | Generate theoretical isotope distributions and cosine similarity scoring |
| [`isotope_pattern_fit_scorer`](tools/metabolomics/spectral_analysis/isotope_pattern_fit_scorer/) | Score isotope pattern fit, detect Cl/Br from M+2 enhancement |
| [`massql_query_tool`](tools/metabolomics/spectral_analysis/massql_query_tool/) | Query mzML data using MassQL-like syntax |

#### Compound Annotation (4 tools)

| Tool | Description |
|------|-------------|
| [`van_krevelen_data_generator`](tools/metabolomics/compound_annotation/van_krevelen_data_generator/) | Compute H:C and O:C ratios, classify into biochemical compound classes |
| [`kendrick_mass_defect_analyzer`](tools/metabolomics/compound_annotation/kendrick_mass_defect_analyzer/) | Compute Kendrick mass defect for homologous series detection (CH2, CF2, etc.) |
| [`suspect_screener`](tools/metabolomics/compound_annotation/suspect_screener/) | Match detected masses against suspect screening lists (CompTox, NORMAN) |
| [`metabolite_class_predictor`](tools/metabolomics/compound_annotation/metabolite_class_predictor/) | Predict compound class from mass defect, element ratios, and RDBE |

#### Drug Metabolism (2 tools)

| Tool | Description |
|------|-------------|
| [`drug_metabolite_screener`](tools/metabolomics/drug_metabolism/drug_metabolite_screener/) | Predict Phase I/II drug metabolites and screen mzML for matches |
| [`mass_difference_network_builder`](tools/metabolomics/drug_metabolism/mass_difference_network_builder/) | Connect features by known biotransformation mass differences |

#### Isotope Labeling (2 tools)

| Tool | Description |
|------|-------------|
| [`isotope_label_detector`](tools/metabolomics/isotope_labeling/isotope_label_detector/) | Detect 13C/15N-labeled metabolites by paired feature analysis |
| [`mid_natural_abundance_corrector`](tools/metabolomics/isotope_labeling/mid_natural_abundance_corrector/) | Correct mass isotopomer distributions for natural 13C abundance |

#### Lipidomics (2 tools)

| Tool | Description |
|------|-------------|
| [`lipid_species_resolver`](tools/metabolomics/lipidomics/lipid_species_resolver/) | Enumerate acyl chain combinations from sum-composition lipid annotations |
| [`lipid_ecn_rt_predictor`](tools/metabolomics/lipidomics/lipid_ecn_rt_predictor/) | Predict lipid RT from Equivalent Carbon Number |

#### Export (3 tools)

| Tool | Description |
|------|-------------|
| [`gnps_fbmn_exporter`](tools/metabolomics/export/gnps_fbmn_exporter/) | Export MS2 + quantification in GNPS Feature-Based Molecular Networking format |
| [`sirius_exporter`](tools/metabolomics/export/sirius_exporter/) | Export features + MS2 data to SIRIUS .ms format |
| [`kovats_ri_calculator`](tools/metabolomics/export/kovats_ri_calculator/) | Calculate Kovats Retention Index from alkane standards for GC-MS |

---

## License

BSD 3-Clause — see [LICENSE](LICENSE).
