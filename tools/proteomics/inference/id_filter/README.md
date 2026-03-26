# ID Filter

Filter peptide/protein identifications by score threshold and decoy status.

## Usage

```bash
# Filter by score threshold
python id_filter.py --input peptides.idXML --output filtered.idXML --score-threshold 0.05

# Remove decoy hits
python id_filter.py --input peptides.idXML --output filtered.idXML --remove-decoys

# Combined filtering
python id_filter.py --input peptides.idXML --output filtered.idXML --score-threshold 0.05 --remove-decoys
```

## Requirements

```
pip install pyopenms click
```

## How it works

Wraps `pyopenms.IDFilter` to filter peptide identifications. Supports score-based filtering (`filterHitsByScore`) and decoy removal (`removeDecoyHits`). Reports before/after hit counts.
