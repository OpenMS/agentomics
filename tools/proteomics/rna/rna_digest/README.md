# RNA Digest

In silico RNA digestion with common RNases.

## Supported Enzymes

- **RNase_T1** - cleaves after G
- **RNase_A** - cleaves after C and U (pyrimidines)
- **RNase_T2** - cleaves after any nucleotide
- **Cusativin** - cleaves after C

## Usage

```bash
python rna_digest.py --sequence AAUGCAAUGG --enzyme RNase_T1
python rna_digest.py --sequence AAUGCAAUGG --enzyme RNase_A --missed-cleavages 1 --output fragments.tsv
```

## Options

- `--sequence` - RNA sequence (A, C, G, U characters)
- `--enzyme` - RNase enzyme name
- `--missed-cleavages` - Maximum missed cleavages (default: 0)
- `--output` - Output TSV file (optional)
