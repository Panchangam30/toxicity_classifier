# Input Preprocessing Module
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import AllChem


def process_input(smiles):
    """
    Convert SMILES to RDKit molecule. Compute MACCS keys and ECFP (Morgan) fingerprints.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        # Invalid SMILES
        return {"molecule": None, "descriptors": {}, "error": "Invalid SMILES string."}

    # MACCS keys (166 bits)
    maccs_fp = MACCSkeys.GenMACCSKeys(mol)
    maccs_bits = list(maccs_fp)

    # ECFP (Morgan fingerprint, radius=2, 2048 bits)
    ecfp_fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    ecfp_bits = list(ecfp_fp)

    descriptors = {
        "MACCS": maccs_bits,
        "ECFP": ecfp_bits
    }
    return {"molecule": mol, "descriptors": descriptors} 