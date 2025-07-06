# hERG Blocker Model

**Model:** hERG_Karim-CNN

## Description

This model predicts the likelihood that a compound blocks the hERG potassium channel, which is a big risk for cardiac toxicity.

## Dependencies

- tdc (and its dependencies)

Install with:

```bash
pip install tdc
```

## How to Run

```bash
python models/run_Pred-hERG.py
```

You will be prompted to enter a SMILES string.

## Input

- **SMILES string** (e.g., `CCO`)

## Output

A JSON file `examples/output_hERG.json` with:

```json
{
  "model_name": "hERG_Karim-CNN",
  "score": 0.51,
  "flag": true
}
```

## Threshold

- **hERG Blocker if:** `score >= 0.5`

## Model Download & Caching

- The model is automatically downloaded and cached in a `.tdc_cache/` directory in your project root.

## Example Usage

```bash
$ python models/run_Pred-hERG.py
🔬 Enter a SMILES string: CCO
{
  "model_name": "hERG_Karim-CNN",
  "score": 0.51,
  "flag": true
}
Result written be in examples/output_hERG.json
```
