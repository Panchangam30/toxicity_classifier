import json, sys, traceback, os
from tdc import tdc_hf_interface

MODEL_ID = "hERG_Karim-CNN"

def predict(smiles: str):
    tdc_hf = tdc_hf_interface(MODEL_ID)
    cache_dir = os.path.join(os.path.dirname(__file__), '..', '.tdc_cache')
    os.makedirs(cache_dir, exist_ok=True)
    dp_model = tdc_hf.load_deeppurpose(cache_dir)
    result = tdc_hf.predict_deeppurpose(dp_model, [smiles])
    return result

def main():
    try:
        smiles = input("Enter a SMILES string: ").strip()
        if not smiles:
            sys.exit("No SMILES string provided.")
        score = predict(smiles)
        print(json.dumps({
            "model_name": MODEL_ID,
            "score": score,
            "flag": score >= 0.5,
            "organ": "heart"
        }, indent=2))
        output_path = "examples/output_hERG.json"
        with open(output_path, "w") as f:
            json.dump({
                "model_name": MODEL_ID,
                "score": score,
                "flag": score >= 0.5,
                "organ": "heart"
            }, f, indent=2)
        print(f"Result written to {output_path}")
    except Exception as e:
        traceback.print_exc()
        sys.exit(f"inference failed: {e}")

if __name__ == "__main__":
    main()
