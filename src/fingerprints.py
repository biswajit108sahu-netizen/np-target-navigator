import pickle

import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def smiles_to_fp(smiles):
    """Morgan fingerprint for a SMILES string, or None if invalid."""
    mol = Chem.MolFromSmiles(smiles)
    return _gen.GetFingerprint(mol) if mol is not None else None


def build_ligand_library(
    in_path="data/processed/activities_clean.parquet",
    out_path="data/processed/ligand_fps.pkl",
):
    df = pd.read_parquet(in_path, columns=["smiles_std"])
    smiles_list = df["smiles_std"].drop_duplicates().tolist()
    fps = {}
    for s in tqdm(smiles_list, desc="Fingerprinting"):
        fp = smiles_to_fp(s)
        if fp is not None:
            fps[s] = fp
    with open(out_path, "wb") as f:
        pickle.dump(fps, f)
    print(f"Saved {len(fps):,} fingerprints")


if __name__ == "__main__":
    build_ligand_library()