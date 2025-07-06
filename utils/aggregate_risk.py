import os
import json
from typing import Dict, Any, List

def get_toxicity_summary(input_folder: str) -> Dict[str, Any]:
    
    results: List[Dict[str, Any]] = []
    flagged_models = 0
    organ_risks = {}
    low_confidence = []

    for fname in os.listdir(input_folder):
        if fname.endswith('.json') and fname not in ('toxicity_summary.json', 'toxicitysummary.json'):
            with open(os.path.join(input_folder, fname), 'r') as f:
                result = json.load(f)
                if not isinstance(result, dict):
                    continue  
                results.append(result)
                if result.get("flag"):
                    flagged_models += 1
                    
                    if "organ" in result:
                        organ = result["organ"]
                        organ_risks.setdefault(organ, []).append(result["model_name"])
                
                score = result.get("score")
                if score is not None and 0.45 < score < 0.55:
                    low_confidence.append(result["model_name"])

    toxicity_flag = flagged_models >= 3

    summary = {
        "model_results": results,
        "flagged_models": flagged_models,
        "toxicity_flag": toxicity_flag,
        "organ_risks": organ_risks,
        "low_confidence_models": low_confidence
    }
    
    out_path = os.path.join(input_folder, "toxicity_summary.json")
    with open(out_path, 'w') as f:
        json.dump(summary, f, indent=2)
    return summary

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python aggregate_risk.py <results_folder>")
        sys.exit(1)
    folder = sys.argv[1]
    summary = get_toxicity_summary(folder)
    print(json.dumps(summary, indent=2)) 