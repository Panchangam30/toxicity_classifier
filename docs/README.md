# Toxicity & Safety Pipeline

## Overview

A modular pipeline to evaluate the toxicity and safety profile of small-molecule drug candidates. The pipeline flags safety concerns, produces a composite toxicity score, and surfaces interpretability outputs (risk breakdown, confidence).

## Features

- Input: SMILES string (required), optional image/3D descriptors for future modules
- Output: JSON with composite score, per-module risks, flags, and confidence
- Modular: Each module is a separate script/function
- Stubs: Where real models are not yet integrated, stubs ensure the pipeline runs end-to-end

## Quickstart

### 1. Clone the repository

```
git clone <your-repo-url>
cd toxicity_classifier
```

### 2. Set up environment (recommended: conda)

```
conda create -n toxenv python=3.9
conda activate toxenv
conda install -c conda-forge rdkit
pip install -r requirements.txt
```

### 3. Model files

Place the following model files in the project root:

- `carcinogen_classifier_rf.joblib`
- `toxcast_classifier_xgb.joblib`
  (Train or request these if not present.)

### 4. Run the pipeline

```
python tox_pipeline/run_pipeline.py
```

Enter a SMILES string when prompted. Example:

```
c1ccc(cc1)[N+](=O)[O-]
```

## Output

A JSON file (`tox_pipeline/toxicity_summary.json`) and console output with fields:

- composite_score
- organ_toxicity, neurotoxicity, mitochondrial_toxicity, tissue_accumulation, morphological_cytotoxicity, immunotoxicity
- structural_alerts, carcinogenicity, toxcast, flags, model_confidence

## Modules

- Input Preprocessing
- Structural Alerts (PAINS/BRENK, SMARTS-based)
- General Toxicity (Carcinogenicity, ToxCast: real models; LD50: TODO)
- Organ, Neuro, Mitochondrial, Tissue Accumulation, Morphological, Immunotoxicity (stubs)
- Explainability & Confidence (stub)
- Scoring & Aggregation

## Notes

- Stubs are clearly marked with TODOs for future model integration
- For best results, use Python 3.8–3.10
- For issues, see module-specific READMEs in `docs/`
