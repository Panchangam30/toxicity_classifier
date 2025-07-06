import json, sys, torch, traceback
from transformers import AutoTokenizer, AutoModelForSequenceClassification

MODEL_ID = "mikemayuare/SMILY-BPE-tox21"

LABELS = [
    "NR-AR","NR-AR-LBD","NR-AhR","NR-Aromatase","NR-ER","NR-ER-LBD",
    "NR-PPAR-γ","SR-ARE","SR-ATAD5","SR-HSE","SR-MMP","SR-p53"
]
THRESHOLD = 0.50

tokenizer = AutoTokenizer.from_pretrained(MODEL_ID)
model     = AutoModelForSequenceClassification.from_pretrained(MODEL_ID)
model.eval()

def predict(smiles: str) -> float:
    with torch.no_grad():
        tokens = tokenizer(smiles, return_tensors="pt", truncation=True, padding=True)
        logits = model(**tokens).logits                      
        probs  = torch.sigmoid(logits).squeeze().tolist()    
        return max(probs)                                    

def main() -> None:
    try:
        smiles = input("Enter a SMILES string: ").strip()
        if not smiles:
            sys.exit("No SMILES string provided.")
        score = round(predict(smiles), 4)
        print(json.dumps({
            "model_name": MODEL_ID,
            "score": score,
            "flag": score >= THRESHOLD,
            "organ": "general"
        }, indent=2))
        output_path = "examples/output_DeepTox.json"
        with open(output_path, "w") as f:
            json.dump({
                "model_name": MODEL_ID,
                "score": score,
                "flag": score >= THRESHOLD,
                "organ": "general"
            }, f, indent=2)
        print(f"Result written to {output_path}")
    except Exception as e:
        traceback.print_exc()
        sys.exit(f"inference failed: {e}")

if __name__ == "__main__":
    main()
