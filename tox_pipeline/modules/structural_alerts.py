# Structural Alerts Module using pure-Python SMARTS-based PAINS and BRENK filters
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# Minimal PAINS and BRENK SMARTS patterns for demonstration
PAINS_SMARTS = [
    ("PAINS_A", "[O;D2]-[C;D3](=O)-[C;D3](=O)-[O;D2]"),
    ("PAINS_B", "[N;D2]=[N;D2]")
    # Add more PAINS patterns as needed
]
BRENK_SMARTS = [
    ("BRENK_A", "[N+](=O)[O-]"),  # Nitro group
    ("BRENK_B", "C(=O)N[NH2]")    # Hydrazide
    # Add more BRENK patterns as needed
]

# Pre-compile SMARTS
PAINS_PATTERNS = [(name, Chem.MolFromSmarts(sma)) for name, sma in PAINS_SMARTS]
BRENK_PATTERNS = [(name, Chem.MolFromSmarts(sma)) for name, sma in BRENK_SMARTS]

def check_structural_alerts(smiles):
    """
    Check for PAINS and BRENK structural alerts using SMARTS patterns (pure Python).
    Returns a list of triggered alert names.
    """
    alerts = []
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print("[WARNING] Invalid SMILES for structural alerts.")
            return []
        # PAINS
        for name, patt in PAINS_PATTERNS:
            if patt is not None and mol.HasSubstructMatch(patt):
                alerts.append(f"PAINS: {name}")
        # BRENK
        for name, patt in BRENK_PATTERNS:
            if patt is not None and mol.HasSubstructMatch(patt):
                alerts.append(f"BRENK: {name}")
        return alerts
    except Exception as e:
        print(f"[WARNING] Structural alerts could not be applied: {e}")
        return [] 