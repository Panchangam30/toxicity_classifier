# Input Preprocessing Module

## Purpose

Converts a SMILES string into an RDKit molecule and computes molecular descriptors for downstream models.

## Features

- Converts SMILES to RDKit Mol object
- Computes MACCS keys (166-bit) and ECFP (Morgan, 2048-bit) fingerprints
- Returns a dictionary with molecule and descriptors
- Handles invalid SMILES gracefully
- 3D conformer generation is not yet implemented (TODO)

## Usage

Called automatically by the main pipeline. Example:

```python
from modules.input_preprocessing import process_input
result = process_input('CCO')
print(result)
```

## Output

- `molecule`: RDKit Mol object
- `descriptors`: dict with 'MACCS' and 'ECFP' keys
- `error`: present if SMILES is invalid
