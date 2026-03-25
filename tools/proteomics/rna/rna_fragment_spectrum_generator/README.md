# RNA Fragment Spectrum Generator

Generate theoretical RNA fragment spectra including c, y, w, and a-B ion series.

## Ion Types

- **c ions** - 5' fragments from 3'-P-O bond cleavage
- **y ions** - 3' complementary fragments
- **w ions** - 3' fragments with base loss
- **a-B ions** - 5' fragments with base loss

## Usage

```bash
python rna_fragment_spectrum_generator.py --sequence AAUGC --charge 2
python rna_fragment_spectrum_generator.py --sequence AAUGC --charge 1 --output fragments.tsv
```

## Options

- `--sequence` - RNA sequence (A, C, G, U)
- `--charge` - Charge state (default: 1)
- `--output` - Output TSV file (optional)
