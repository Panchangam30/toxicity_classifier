# Toxicity & Safety Pipeline

## Project Overview

This project provides a modular, extensible pipeline for evaluating the toxicity and safety profile of small-molecule drug candidates. It is designed to help researchers, computational chemists, and drug developers flag safety concerns early, prioritize compounds, and gain interpretable insights into toxicity risks.

---

## What Problems Does This Project Solve?

- **Early Risk Assessment:** Flags potential toxicity and safety issues before costly experiments or clinical trials.
- **Composite Scoring:** Aggregates multiple toxicity endpoints into a single, interpretable score.
- **Interpretability:** Surfaces risk breakdowns, confidence scores, and high-risk flags for actionable insights.
- **Modularity:** Each toxicity endpoint is a separate module, making it easy to swap in new models or add new endpoints.
- **Reproducibility:** Provides a clear, documented workflow for consistent toxicity evaluation.

---

## Who Is This For?

- Drug discovery scientists
- Computational chemists
- Toxicologists
- Data scientists in pharma/biotech
- Anyone needing rapid, interpretable toxicity profiling of small molecules

---

## Scientific Context

Toxicity prediction is a critical step in drug discovery. This pipeline brings together state-of-the-art models and best practices from cheminformatics, machine learning, and toxicology, including:

- **RDKit** for molecular representation and descriptor calculation.
- **Rule-based filters** (e.g., PAINS, BRENK) for structural alerts.
- **Machine learning models** for endpoints like carcinogenicity and ToxCast (real models integrated; others are stubbed).
- **Stubs/mocks** for rapid prototyping and pipeline development, with clear TODOs for future model integration.

---

## What Does the Pipeline Do?

- **Input:** Accepts a SMILES string (chemical structure) as input.
- **Processing:**
  - Converts SMILES to a molecule object using RDKit.
  - Computes molecular descriptors (MACCS, ECFP).
- **Module Evaluation:**
  - Runs a series of modules (general toxicity, toxcast, organ-specific, neurotoxicity, etc.).
  - **General toxicity and toxcast modules use real trained models.**
  - Other modules are currently stubbed/mocked (return random or fixed values).
- **Aggregation:**
  - Combines all module outputs into a composite toxicity score.
  - Adds interpretability outputs (flags, confidence, risk breakdown).
- **Output:**
  - Prints a complete JSON report to the terminal.
  - Saves the report to `tox_pipeline/toxicity_summary.json`.

---

## Models and Module Types

- **Input Preprocessing:** Uses RDKit to convert SMILES to molecules and compute descriptors (MACCS, ECFP).
- **Structural Alerts:** (Stubbed) Would use RDKit PAINS/BRENK filters; currently returns no alerts due to environment limitations.
- **General Toxicity:** Uses a real trained model (Random Forest, XGBoost, or LightGBM) for carcinogenicity (Carcinogens_Lagunin dataset).
- **ToxCast:** Uses a real trained model (XGBoost, LightGBM, or Random Forest) for a selected ToxCast endpoint.
- **Organ-Specific Toxicity:** Placeholder for ML models (e.g., H-optimus-0, UNI, Merlin); currently mocked.
- **Neurotoxicity, Mitochondrial Toxicity, Morphological Cytotoxicity, Immunotoxicity:** All are modular and ready for real models; currently mocked.
- **Tissue Accumulation:** Placeholder for pkCSM or similar models; currently mocked.
- **Explainability & Confidence:** Placeholder for ensemble or calibration-based confidence; currently mocked.
- **Scoring & Aggregation:** Implements a weighted formula to combine module outputs and flag high-risk endpoints.

**Note:** Only the general toxicity and toxcast modules currently use real models. All other modules are stubbed/mocked for demonstration. Real models can be integrated as needed.

---

## Model Files

- Trained model files (e.g., `carcinogen_classifier_rf.joblib`, `toxcast_classifier_xgb.joblib`) are saved in the project root after running the training script (`general_toxicity_train.py`).
- The pipeline expects these files to be present in the project root. If missing, run:
  ```sh
  python general_toxicity_train.py
  ```
- You can swap in new models by updating the filenames in `tox_pipeline/modules/general_toxicity.py`.

---

## Example Workflow

1. Researcher provides a SMILES string for a new compound.
2. Pipeline processes the input, computes descriptors, and runs all toxicity modules.
3. **General toxicity and toxcast modules use real trained models to make predictions.**
4. Each other module returns a risk score (mocked or real, as available).
5. The scoring module aggregates results into a composite score and flags high-risk endpoints.
6. The output JSON provides a full toxicity profile, ready for review or further analysis.

---

## How to Run the Pipeline

1. **Set up your environment:**
   - Recommended: Use conda
   - Create and activate environment:
     ```sh
     conda create -n toxenv python=3.10 rdkit numpy -c conda-forge
     conda activate toxenv
     ```
2. **Install dependencies (if not using conda):**
   ```sh
   pip install -r requirements.txt
   ```
3. **Train the models (if not already present):**
   ```sh
   python general_toxicity_train.py
   ```
4. **Run the pipeline:**
   ```sh
   python tox_pipeline/run_pipeline.py
   ```
   - The output will be printed to the terminal and saved to `tox_pipeline/toxicity_summary.json`.

---

## Output Fields

- `composite_score`: Overall toxicity risk (0–1)
- `organ_toxicity`: Per-organ risk scores (stubbed)
- `neurotoxicity`, `mitochondrial_toxicity`, `morphological_cytotoxicity`, `immunotoxicity`: Module risk scores (stubbed)
- `tissue_accumulation`: Predicted accumulation in liver/brain (stubbed)
- `structural_alerts`: List of triggered alerts (stubbed)
- `carcinogenicity`: Real model prediction and probability
- `toxcast`: Real model prediction and probability
- `flags`: High-risk warnings
- `model_confidence`: Confidence score (stubbed)

---

## Project Structure

- `tox_pipeline/` — Main pipeline package
  - `modules/` — Each toxicity/safety module as its own script (input preprocessing, structural alerts, general toxicity, etc.)
  - `utils/` — Helper functions (e.g., error handling, formatting)
  - `run_pipeline.py` — Main entry script to run the pipeline
  - `toxicity_summary.json` — Output file with the latest pipeline results
- `docs/` — Documentation for each module and pipeline usage
  - `README.md` — Documentation index
  - `README_run_pipeline.md` — How to run the pipeline
  - `README_<module>.md` — Details for each module (input preprocessing, general toxicity, etc.)
- `requirements.txt` — Python dependencies for running the pipeline
- `examples/` — Example input and output files for testing and demonstration
- `models/` — (If present) Scripts for training or running real models
- `tests/` — Test scripts for validating pipeline functionality

---

## Extending the Pipeline

- **Swap in real models:** Replace stubs with trained models for any module.
- **Add new endpoints:** Create new modules for additional toxicity endpoints or data types.
- **Integrate with other tools:** Use the output JSON for downstream analysis, visualization, or reporting.
- **Collaborate:** The modular structure and documentation make it easy for teams to contribute.

---

## Future Directions

- Integrate real machine learning models for each endpoint
- Add support for image-based, omics, or 3D structure data
- Improve explainability and confidence estimation
- Add a web or CLI interface for batch processing
- Expand documentation with example outputs and troubleshooting

---

## Where to Find More Information

- See the `docs/` folder for detailed documentation on each module and how to run the pipeline.
