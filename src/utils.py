from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

_uncharger = rdMolStandardize.Uncharger()

def standardize(smiles):
    """Return a cleaned canonical SMILES, or None if invalid."""
    if not isinstance(smiles, str):
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol = rdMolStandardize.Cleanup(mol)
        mol = rdMolStandardize.FragmentParent(mol)
        mol = _uncharger.uncharge(mol)
        return Chem.MolToSmiles(mol)
    except Exception:
        return None