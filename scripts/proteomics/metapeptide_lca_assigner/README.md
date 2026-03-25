# Metapeptide LCA Assigner

Compute lowest common ancestor taxonomy from peptide-protein mappings.

## Usage

```bash
python metapeptide_lca_assigner.py --peptides peptides.tsv --fasta metadb.fasta --taxonomy lineage.tsv --output taxonomy.tsv
```

## Input Format

- `peptides.tsv`: column `peptide`
- `metadb.fasta`: Meta-proteomics FASTA database
- `lineage.tsv`: columns `protein`, `lineage` (semicolon-separated taxonomy)

## Output

- `taxonomy.tsv` - Per-peptide LCA assignments with depth and specificity
