import json
from typing import Dict, Any, List

def format_output(final_summary: Dict[str, Any]) -> Dict[str, Any]:

    toxicity_flag = final_summary.get("toxicity_flag", False)
    organ_systems_flagged = list(final_summary.get("organ_risks", {}).keys())
    num_models_flagged = final_summary.get("flagged_models", 0)
    
    model_scores = {
        r["model_name"]: r["score"]
        for r in final_summary.get("model_results", [])
        if isinstance(r, dict) and "model_name" in r and "score" in r
    }
    
    scores = list(model_scores.values())
    n_flagged = num_models_flagged
    n_models = len(model_scores)
    
    low_conf = False
    if n_models > 1:
        if 0 < n_flagged < n_models:  
            low_conf = True
        near_half = sum(1 for s in scores if isinstance(s, (float, int)) and 0.45 < s < 0.55)
        if near_half >= n_models // 2:
            low_conf = True
    result = {
        "toxicity_flag": toxicity_flag,
        "organ_systems_flagged": organ_systems_flagged,
        "num_models_flagged": num_models_flagged,
        "model_scores": model_scores
    }
    if low_conf:
        result["confidence_flag"] = "low"
    return result

def validate_schema(json_dict: Dict[str, Any]) -> None:
    """
    Validates the toxicity report schema. Raises ValueError if invalid.
    """
    if not isinstance(json_dict, dict):
        raise ValueError("Report must be a dictionary.")
    if not isinstance(json_dict.get("toxicity_flag"), bool):
        raise ValueError("toxicity_flag must be a boolean.")
    if not isinstance(json_dict.get("organ_systems_flagged"), list):
        raise ValueError("organ_systems_flagged must be a list.")
    if not isinstance(json_dict.get("num_models_flagged"), int):
        raise ValueError("num_models_flagged must be an integer.")
    if not isinstance(json_dict.get("model_scores"), dict):
        raise ValueError("model_scores must be a dictionary.")
    for k, v in json_dict["model_scores"].items():
        if not isinstance(k, str):
            raise ValueError("All model names in model_scores must be strings.")
        if v is not None and not (isinstance(v, (float, int, bool))):
            raise ValueError("All model scores must be float, int, bool, or None.")
