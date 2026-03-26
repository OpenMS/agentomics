# Decoy Database Generator

Generate a target-decoy FASTA database by reversing or shuffling protein sequences. The output contains all original (target) entries followed by the generated decoy entries with a `DECOY_` identifier prefix.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

### Command Line

```bash
# Reverse decoys (default)
python decoy_database_generator.py --input target.fasta --output target_decoy.fasta --method reverse

# Shuffled decoys
python decoy_database_generator.py --input target.fasta --output target_decoy.fasta --method shuffle
```

### Python API

```python
from decoy_database_generator import generate_decoys

count = generate_decoys("target.fasta", "target_decoy.fasta", method="reverse")
print(f"Generated {count} decoy proteins")
```

## Options

| Option     | Default   | Description                          |
|------------|-----------|--------------------------------------|
| `--input`  | required  | Input target FASTA file              |
| `--output` | required  | Output target-decoy FASTA file       |
| `--method` | `reverse` | Decoy method: `reverse` or `shuffle` |

## Output

A FASTA file containing the original target proteins followed by decoy proteins. Each decoy identifier is prefixed with `DECOY_`.
