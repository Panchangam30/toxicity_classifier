# General Toxicity Module
import joblib
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import os

# --- Featurization ---
def smiles_to_ecfp(smiles, radius=2, n_bits=2048):
    """
    Convert a SMILES string to an ECFP (Morgan) fingerprint.
    Returns a numpy array of length n_bits.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(n_bits)
    return np.array(AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits))

# --- Model Loading ---
# Always resolve model paths from the project root
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../'))
CARCINOGEN_MODEL_PATH = os.path.join(PROJECT_ROOT, 'carcinogen_classifier_rf.joblib')
TOXCAST_MODEL_PATH = os.path.join(PROJECT_ROOT, 'toxcast_classifier_xgb.joblib')

carc_model = joblib.load(CARCINOGEN_MODEL_PATH)
toxcast_model = joblib.load(TOXCAST_MODEL_PATH)

# TODO: Add model loading for additional general toxicity endpoints here (e.g., LD50)
# Example:
# LD50_MODEL_PATH = os.path.join(PROJECT_ROOT, 'ld50_model.joblib')
# ld50_model = joblib.load(LD50_MODEL_PATH)

# --- Prediction Functions ---
def predict_carcinogenicity(smiles: str):
    """
    Predict carcinogenicity for a given SMILES string.
    Returns a dict with prediction (class) and probability (confidence).
    """
    features = smiles_to_ecfp(smiles)
    pred = carc_model.predict([features])[0]
    prob = carc_model.predict_proba([features])[0]
    return {'prediction': int(pred), 'probability': float(np.max(prob))}


def predict_toxcast(smiles: str):
    """
    Predict ToxCast toxicity for a given SMILES string.
    Returns a dict with prediction (class) and probability (confidence).
    """
    features = smiles_to_ecfp(smiles)
    pred = toxcast_model.predict([features])[0]
    prob = toxcast_model.predict_proba([features])[0]
    return {'prediction': int(pred), 'probability': float(np.max(prob))}

# TODO: Add prediction functions for additional endpoints here (e.g., LD50)
# def predict_ld50(smiles: str):
#     """
#     Predict LD50 for a given SMILES string using a real trained model.
#     Returns a float value (mg/kg).
#     """
#     features = smiles_to_ecfp(smiles)
#     ld50 = ld50_model.predict([features])[0]
#     return float(ld50)

# --- Example Usage ---
if __name__ == "__main__":
    test_smiles = "CCO"  # Example: ethanol
    print("Carcinogenicity:", predict_carcinogenicity(test_smiles))
    print("ToxCast:", predict_toxcast(test_smiles))
    # TODO: Add example usage for additional endpoints as implemented 