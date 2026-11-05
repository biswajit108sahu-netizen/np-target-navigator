import pickle
import sys
from functools import lru_cache

import pandas as pd
from rdkit import DataStructs

from fingerprints import smiles_to_fp
from utils import standardize

ACTS_PATH = "data/processed/activities_clean.parquet"
FPS_PATH = "data/processed/ligand_fps.pkl"

CURCUMIN = "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O"


@lru_cache(maxsize=1)
def _load():
    with open(FPS_PATH, "rb") as f:
        fps = pickle.load(f)
    acts = pd.read_parquet(
        ACTS_PATH,
        columns=["smiles_std", "target_chembl_id", "target_name", "pchembl_value"],
    )
    acts["target_name"] = acts["target_name"].fillna("")
    return list(fps.keys()), list(fps.values()), acts


def predict_targets(smiles, top_n=10, k=50, min_sim=0.3, exclude_exact=False):
    """Rank likely targets for a query compound by similarity to known ligands.

    k             number of nearest known ligands to use
    min_sim       ignore neighbours less similar than this (Tanimoto)
    exclude_exact drop near-identical library hits (use for validation)
    The result's .attrs["max_sim"] is the best similarity to the library
    (low values mean the compound is outside the model's domain).
    """
    std = standardize(smiles)
    if std is None:
        raise ValueError("Invalid SMILES")
    qfp = smiles_to_fp(std)

    lib_smiles, lib_fps, acts = _load()
    sims = pd.Series(DataStructs.BulkTanimotoSimilarity(qfp, lib_fps), index=lib_smiles)
    if exclude_exact:
        sims = sims[sims < 0.999]

    max_sim = float(sims.max())
    neigh = sims.nlargest(k)
    neigh = neigh[neigh >= min_sim]

    cols = ["target_chembl_id", "target_name", "score",
            "confidence", "max_sim", "n_neighbors", "mean_pchembl"]
    if neigh.empty:
        out = pd.DataFrame(columns=cols)
        out.attrs["max_sim"] = max_sim
        return out

    hits = acts[acts["smiles_std"].isin(neigh.index)].copy()
    hits["sim"] = hits["smiles_std"].map(neigh)

    out = (
        hits.groupby(["target_chembl_id", "target_name"])
        .agg(score=("sim", "sum"),
             max_sim=("sim", "max"),
             n_neighbors=("smiles_std", "nunique"),
             mean_pchembl=("pchembl_value", "mean"))
        .reset_index()
    )
    out["confidence"] = (out["score"] / neigh.sum()).clip(upper=1.0)
    out = out.sort_values(["score", "max_sim"], ascending=False).head(top_n)
    out = out[cols].reset_index(drop=True)
    out.attrs["max_sim"] = max_sim
    return out


if __name__ == "__main__":
    smi = sys.argv[1] if len(sys.argv) > 1 else CURCUMIN
    res = predict_targets(smi)
    print(f"Best similarity to library: {res.attrs['max_sim']:.2f}")
    print(res.round(2).to_string(index=False))