# Simple Search Engine

Run a simple peptide identification search against an mzML file and a FASTA protein database using the OpenMS `SimpleSearchEngineAlgorithm`. Results are written in idXML format.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

### Command Line

```bash
python simple_search_engine.py --input run.mzML --database proteins.fasta --output results.idXML --precursor-tol 10.0 --fragment-tol 20.0
```

### Python API

```python
from simple_search_engine import search

psm_count = search("run.mzML", "proteins.fasta", "results.idXML", precursor_tol=10.0, fragment_tol=20.0)
print(f"Found {psm_count} PSMs")
```

## Options

| Option            | Default | Description                           |
|-------------------|---------|---------------------------------------|
| `--input`         | required| Input mzML file                       |
| `--database`      | required| FASTA protein database                |
| `--output`        | required| Output idXML file                     |
| `--precursor-tol` | `10.0`  | Precursor mass tolerance in ppm       |
| `--fragment-tol`  | `20.0`  | Fragment mass tolerance in ppm        |

## Output

An idXML file containing protein and peptide identifications from the search.
