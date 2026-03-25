# Search Result Merger

Merge multiple identification TSV files with union or intersection consensus.

## Usage

```bash
python search_result_merger.py --inputs engine1.tsv engine2.tsv --method union --output merged.tsv
python search_result_merger.py --inputs engine1.tsv engine2.tsv --method intersection --output merged.tsv
```

## Methods

- **union** - Include all PSMs from any search engine
- **intersection** - Only include PSMs found in all search engines
