# FASTA Decoy Validator

Check if a FASTA database contains decoy sequences and validate prefix consistency.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# Validate with default DECOY_ prefix
python fasta_decoy_validator.py --input db.fasta

# Validate with custom prefix
python fasta_decoy_validator.py --input db.fasta --decoy-prefix REV_ --output validation.json
```
