from rdkit import Chem
from PIL import Image
import numpy as np

def validate_smiles(smiles):
    """
    Validate a SMILES string. Return canonical SMILES if valid, else None.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)

def preprocess_histology_image(image_path, size=(224, 224), mean=(0.707223, 0.578729, 0.703617), std=(0.211883, 0.230117, 0.177517)):
    """
    Load and preprocess a histology image for H-optimus-0.
    Returns a normalized numpy array.
    """
    img = Image.open(image_path).convert("RGB").resize(size)
    arr = np.array(img).astype(np.float32) / 255.0
    arr = (arr - np.array(mean)) / np.array(std)
    return arr 