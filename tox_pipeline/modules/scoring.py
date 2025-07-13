# Scoring & Aggregation Module

def aggregate_scores(module_outputs):
    """
    Aggregate module outputs into a composite toxicity score (0-1).
    Add flags for high-risk values.
    TODO: Refine normalization and thresholding as needed.
    """
    
    general_tox = module_outputs.get("general_tox", {})
    # Use the 'prediction' value from both carcinogenicity and toxcast, averaged
    carcinogenicity_pred = general_tox.get("carcinogenicity", {}).get("prediction", 0)
    toxcast_pred = general_tox.get("toxcast", {}).get("prediction", 0)
    general_tox_avg = (carcinogenicity_pred + toxcast_pred) / 2

    organ_tox = module_outputs.get("organ_tox", {})
    organ_tox_values = [organ_tox.get("cardiotoxicity", 0), organ_tox.get("hepatotoxicity", 0), organ_tox.get("nephrotoxicity", 0)]
    organ_tox_avg = sum(organ_tox_values) / len(organ_tox_values) if organ_tox_values else 0
    neurotox = module_outputs.get("neurotox", {}).get("neurotoxicity", 0)
    mito_tox = module_outputs.get("mito_tox", {}).get("mitochondrial_toxicity", 0)
    morpho_tox = module_outputs.get("morpho_tox", {}).get("morphological_cytotoxicity", 0)
    accumulation = module_outputs.get("accumulation", {})
    immunotox = module_outputs.get("immunotoxicity", {}).get("immunotoxicity", 0)
    alerts = module_outputs.get("alerts", [])

    accumulation_penalty = 1 if (accumulation.get("liver") == "high" or accumulation.get("brain") == "high") else 0
    structural_alert_penalty = 1 if len(alerts) > 0 else 0

    # Composite score formula
    score = (
        0.15 * general_tox_avg +
        0.2 * organ_tox_avg +
        0.15 * neurotox +
        0.1 * mito_tox +
        0.1 * morpho_tox +
        0.1 * accumulation_penalty +
        0.1 * immunotox +
        0.1 * structural_alert_penalty
    )
    
    score = max(0, min(1, score))

    flags = []
    if organ_tox.get("hepatotoxicity", 0) > 0.8:
        flags.append("high hepatotoxicity")
    if organ_tox.get("cardiotoxicity", 0) > 0.8:
        flags.append("high cardiotoxicity")
    if morpho_tox > 0.8:
        flags.append("morphological concern")
    if accumulation_penalty:
        flags.append("high tissue accumulation")
    if structural_alert_penalty:
        flags.append("structural alert present")
    if immunotox > 0.8:
        flags.append("high immunotoxicity")
    if mito_tox > 0.8:
        flags.append("high mitochondrial toxicity")
    if neurotox > 0.8:
        flags.append("high neurotoxicity")

    return {"composite_score": score, "flags": flags} 