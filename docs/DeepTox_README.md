# DeepTox Model

**Model:** mikemayuare/SMILY-BPE-tox21

## Description

DeepTox is for predicting general toxicity of drug-like compounds. It outputs a probability score indicating the likelihood of toxicity.

## Dependencies

- torch
- transformers

Install with:

```bash
pip install torch transformers
```

## How to Run

```bash
python models/run_DeepTox.py
```

You will be prompted to enter a SMILES string.

## Input

- **SMILES string** (e.g., `CCO`)

## Output

A JSON file `examples/output_DeepTox.json` with:

```json
{
  "model_name": "DeepTox",
  "score": 0.82,
  "flag": true
}
```

## Threshold

- **Toxicity if:** `score >= 0.5`

## Example Usage

```bash
$ python models/run_DeepTox.py
  Enter a SMILES string: CCO
{
  "model_name": "DeepTox",
  "score": 0.82,
  "flag": true
}
Result will be in examples/output_DeepTox.json
```
