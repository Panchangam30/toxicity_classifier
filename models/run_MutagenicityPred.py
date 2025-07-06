import json, sys, traceback
from transformers import pipeline

MODEL_ID = "ML4chemistry/Toxicity_Prediction_of_Ames_test"
THRESHOLD = 0.5


pipe = pipeline("text-classification", model=MODEL_ID)

def predict(smiles: str) -> float:
    result = pipe(smiles)[0]
    
    score = result["score"] if result["label"] == "LABEL_1" else 1 - result["score"]
    return score

def main():
    try:
        smiles = input("Enter a SMILES string: ").strip()
        if not smiles:
            sys.exit("No SMILES string provided.")
        score = round(predict(smiles), 4)
        print(json.dumps({
            "model_name": MODEL_ID,
            "score": score,
            "flag": score >= THRESHOLD,
            "organ": "mutagenicity"
        }, indent=2))
        output_path = "examples/output_Mutagenicity.json"
        with open(output_path, "w") as f:
            json.dump({
                "model_name": MODEL_ID,
                "score": score,
                "flag": score >= THRESHOLD,
                "organ": "mutagenicity"
            }, f, indent=2)
        print(f"Result written to {output_path}")
    except Exception as e:
        traceback.print_exc()
        sys.exit(f"inference failed: {e}")

if __name__ == "__main__":
    main()
