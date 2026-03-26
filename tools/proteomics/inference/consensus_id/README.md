# Consensus ID

Merge peptide identifications from multiple search engines using consensus scoring algorithms.

## Usage

```bash
# Use 'best' score from any engine
python consensus_id.py --inputs search1.idXML --inputs search2.idXML --output consensus.idXML --algorithm best

# Use average scores
python consensus_id.py --inputs search1.idXML --inputs search2.idXML --output consensus.idXML --algorithm average

# Use rank-based scoring
python consensus_id.py --inputs search1.idXML --inputs search2.idXML --output consensus.idXML --algorithm ranks
```

## Requirements

```
pip install pyopenms click
```

## How it works

Wraps `pyopenms.ConsensusIDAlgorithmBest`, `ConsensusIDAlgorithmAverage`, and `ConsensusIDAlgorithmRanks` to combine identifications from multiple search engines into consensus results. Supports three scoring strategies: best (highest confidence from any engine), average (mean across engines), and ranks (rank-based aggregation).
