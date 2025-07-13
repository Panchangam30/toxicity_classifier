# Main Pipeline Script (`run_pipeline.py`)

## Purpose

Runs the full toxicity & safety pipeline on a given SMILES string.

## How to Run

1. Ensure all dependencies are installed (see main README)
2. Ensure required model files are present in the project root
3. Run:
   ```
   python tox_pipeline/run_pipeline.py
   ```
4. Enter a SMILES string when prompted

## Input

- SMILES string (e.g., `c1ccc(cc1)[N+](=O)[O-]`)

## Output

- Prints a JSON summary to the console
- Saves output to `tox_pipeline/toxicity_summary.json`

## Troubleshooting

- **RDKit install errors:** Use conda for best results
- **Missing model files:** Ensure `carcinogen_classifier_rf.joblib` and `toxcast_classifier_xgb.joblib` are in the project root
- **Python version:** Use Python 3.8–3.10 for best compatibility

## See Also

- Module-specific READMEs in `docs/`
- Main README for setup and usage
