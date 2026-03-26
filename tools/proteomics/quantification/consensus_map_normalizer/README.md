# Consensus Map Normalizer

Normalize intensity values across channels in a consensusXML file using median or quantile normalization.

## Usage

```bash
# Median normalization (default)
python consensus_map_normalizer.py --input consensus.consensusXML --output normalized.consensusXML --method median

# Quantile normalization
python consensus_map_normalizer.py --input consensus.consensusXML --output normalized.consensusXML --method quantile
```

## Dependencies

```
pyopenms
click
```
