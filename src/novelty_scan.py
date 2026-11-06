import pandas as pd
from rdkit import DataStructs
from tqdm import tqdm

from fingerprints import smiles_to_fp
from predict import _load


def main(n=400, seed=7):
    _, lib_fps, _ = _load()
    nps = pd.read_parquet("data/processed/np_clean.parquet")
    sample = nps.sample(n=min(n, len(nps)), random_state=seed).copy()

    best = []
    for smi in tqdm(sample["smiles_std"], desc="Scanning COCONUT sample"):
        fp = smiles_to_fp(smi)
        if fp is None:
            best.append(None)
        else:
            best.append(max(DataStructs.BulkTanimotoSimilarity(fp, lib_fps)))
    sample["max_sim_to_chembl"] = best
    sample = sample.dropna(subset=["max_sim_to_chembl"])
    sample = sample.sort_values("max_sim_to_chembl")
    sample.to_csv("docs/coconut_novelty_scan.csv", index=False)

    s = sample["max_sim_to_chembl"]
    print(f"\nCompounds scanned: {len(s)}")
    print(f"Median similarity to nearest ChEMBL ligand: {s.median():.2f}")
    for t in (0.3, 0.4, 0.7):
        print(f"Share below {t}: {(s < t).mean():.0%}")

    cols = [c for c in sample.columns if c not in ("smiles_std",)][:2]
    print("\nCandidates for the demo (similarity 0.30-0.40):")
    mid = sample[(s >= 0.30) & (s < 0.40)].head(8)
    print(mid[cols + ["smiles_std", "max_sim_to_chembl"]].to_string(index=False))
    print("\nCandidates with no usable neighbours (below 0.30):")
    low = sample[s < 0.30].head(3)
    print(low[cols + ["smiles_std", "max_sim_to_chembl"]].to_string(index=False))


if __name__ == "__main__":
    main()