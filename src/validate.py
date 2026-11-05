import random

import pandas as pd
from tqdm import tqdm

from predict import predict_targets, _load


def main(n=300, seed=42, top_n=10):
    _, _, acts = _load()
    truth = acts.groupby("smiles_std")["target_chembl_id"].agg(set)
    popular = set(acts["target_chembl_id"].value_counts().head(top_n).index)

    random.seed(seed)
    queries = random.sample(list(truth.index), n)

    rows = []
    for smi in tqdm(queries, desc="Validating"):
        true_t = truth[smi]
        res = predict_targets(smi, top_n=top_n, exclude_exact=True)
        preds = list(res["target_chembl_id"])
        rows.append({
            "smiles": smi,
            "max_sim": res.attrs["max_sim"],
            "top1": bool(preds[:1] and preds[0] in true_t),
            "top5": bool(set(preds[:5]) & true_t),
            "top10": bool(set(preds[:10]) & true_t),
            "baseline_top10": bool(popular & true_t),
        })

    df = pd.DataFrame(rows)
    df.to_csv("data/processed/validation_results.csv", index=False)

    print(f"\nOverall (n={n}):")
    print(df[["top1", "top5", "top10", "baseline_top10"]].mean().round(3).to_string())

    df["similarity_bin"] = pd.cut(df["max_sim"], [0, 0.4, 0.7, 1.01],
                                  labels=["<0.4", "0.4-0.7", ">=0.7"])
    print("\nBy similarity to nearest other library compound:")
    print(df.groupby("similarity_bin", observed=True)[
        ["top1", "top5", "top10", "baseline_top10"]].mean().round(3).to_string())
    print(df["similarity_bin"].value_counts().sort_index().to_string())


if __name__ == "__main__":
    main()