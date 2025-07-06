# Mutagenicity (Ames) Model

**Model:** ML4chemistry/Toxicity_Prediction_of_Ames_test

## Description

This model predicts the mutagenicity (Ames test) of drug-like compounds. It outputs a probability score for mutagenicity.

## Dependencies

- transformers

Install with:

```bash
pip install transformers
```

## How to Run

```bash
python models/run_MutagenicityPred.py
```

You will be prompted to enter a SMILES string.

## Input

- **SMILES string** (e.g., `CCO`)

## Output

A JSON file `examples/output_Mutagenicity.json` with:

```json
{
  "model_name": "ML4chemistry/Toxicity_Prediction_of_Ames_test",
  "score": 0.67,
  "flag": true
}
```

## Threshold

- **Mutagenic if:** `score >= 0.5`

## Example Usage

```bash
$ python models/run_MutagenicityPred.py
  Enter a SMILES string: CCO
{
  "model_name": "ML4chemistry/Toxicity_Prediction_of_Ames_test",
  "score": 0.67,
  "flag": true
}
Result written be in examples/output_Mutagenicity.json
```
