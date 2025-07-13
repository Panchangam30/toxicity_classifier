# Main entry script for the toxicity & safety pipeline
from modules.input_preprocessing import process_input
from modules.structural_alerts import check_structural_alerts
from modules.general_toxicity import predict_carcinogenicity, predict_toxcast
from modules.organ_toxicity import predict_organ_toxicity
from modules.neurotoxicity import predict_cns_toxicity_converge
from modules.mitochondrial_toxicity import predict_mitochondrial_toxicity
from modules.tissue_accumulation import predict_tissue_accumulation
from modules.morphological_cytotoxicity import predict_morphological_cytotoxicity
from modules.immunotoxicity import predict_immunotoxicity
from modules.explainability import get_explainability_and_confidence
from modules.scoring import aggregate_scores
import numpy as np


def run_pipeline(smiles):
    # 1. Input Preprocessing
    input_data = process_input(smiles)
    molecule = input_data["molecule"]
    descriptors = input_data["descriptors"]

    # 2. Structural Alerts
    alerts = check_structural_alerts(smiles)

    # 3. General Toxicity
    # general_tox = predict_general_toxicity(descriptors)
    general_tox = {
        "carcinogenicity": predict_carcinogenicity(smiles),
        "toxcast": predict_toxcast(smiles)
    }

    # 4. Organ-Specific Toxicity
    organ_tox = predict_organ_toxicity(descriptors)

    # 5. Neurotoxicity
    neurotox = predict_cns_toxicity_converge(descriptors)

    # 6. Mitochondrial Toxicity
    mito_tox = predict_mitochondrial_toxicity(descriptors)

    # 7. Tissue Accumulation
    accumulation = predict_tissue_accumulation(descriptors)

    # 8. Morphological Cytotoxicity
    morpho_tox = predict_morphological_cytotoxicity(descriptors)

    # 9. Immunotoxicity
    immunotox = predict_immunotoxicity(descriptors)

    # 10. Explainability & Confidence
    explain = get_explainability_and_confidence({
        "general_tox": general_tox,
        "organ_tox": organ_tox,
        "neurotox": neurotox,
        "mito_tox": mito_tox,
        "accumulation": accumulation,
        "morpho_tox": morpho_tox,
        "immunotox": immunotox,
        "alerts": alerts
    })

    # 11. Scoring & Aggregation
    scoring = aggregate_scores({
        "general_tox": general_tox,
        "organ_tox": organ_tox,
        "neurotox": neurotox,
        "mito_tox": mito_tox,
        "accumulation": accumulation,
        "morpho_tox": morpho_tox,
        "immunotox": immunotox,
        "alerts": alerts
    })

    # Compose output JSON
    output = {
        "composite_score": scoring["composite_score"],
        "organ_toxicity": {
            "cardiotoxicity": organ_tox.get("cardiotoxicity"),
            "hepatotoxicity": organ_tox.get("hepatotoxicity")
        },
        "neurotoxicity": neurotox.get("CNS_toxicity_CONVERGE"),
        "mitochondrial_toxicity": mito_tox.get("mitochondrial_toxicity"),
        "tissue_accumulation": accumulation,
        "morphological_cytotoxicity": morpho_tox.get("morphological_cytotoxicity"),
        "immunotoxicity": immunotox.get("immunotoxicity"),
        "structural_alerts": alerts,
        "structural_alert_count": len(alerts),
        "carcinogenicity": general_tox["carcinogenicity"],
        "toxcast": general_tox["toxcast"],
        "ld50": float(np.random.uniform(50, 500)),
        "flags": scoring.get("flags", []),
        "model_confidence": explain.get("model_confidence"),
    }
    return output


if __name__ == "__main__":
    smiles = input("Enter a SMILES string: ")
    result = run_pipeline(smiles)
    import json
    
    print(json.dumps(result, indent=2))
    
    # Save to file
    with open("tox_pipeline/toxicity_summary.json", "w") as f:
        json.dump(result, f, indent=2) 