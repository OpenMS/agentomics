# DIA Scorer

Score DIA spectra against transition lists using DIAScoring isotope-based methods.

## Usage

```bash
python dia_scorer.py --input dia.mzML --transitions transitions.tsv --output scores.tsv
```

## Output

Tab-separated file with columns:
- `transition_group_id` - Transition group identifier
- `peptide_sequence` - Peptide sequence
- `protein_name` - Protein name
- `precursor_mz` - Precursor m/z
- `n_transitions` - Number of transitions in group
- `best_rt` - RT of best-scoring spectrum
- `dotprod_score` - Dot product score
- `manhattan_score` - Manhattan distance score
- `scored` - Whether scoring succeeded

## Dependencies

- pyopenms
- click
