# Sequence Tag Generator

Generate de novo sequence tags from MS2 spectra by matching peak mass differences to amino acid residue masses.

## Usage

```bash
python sequence_tag_generator.py --mz-list "200.1,313.2,426.3,539.4" --intensities "100,200,150,300" \
    --tolerance 0.02 --min-tag-length 3 --output tags.tsv
```
