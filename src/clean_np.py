import pandas as pd
from tqdm import tqdm

from utils import standardize

IN_PATH = "data/raw/coconut.csv"
OUT_PATH = "data/processed/np_clean.parquet"
KEEP_HINTS = ("identifier", "coconut_id", "name", "organism", "formula")


def main():
    with open(IN_PATH, encoding="utf-8", errors="replace") as f:
        header = f.readline()
    sep = "\t" if header.count("\t") > header.count(",") else ","

    cols = list(pd.read_csv(IN_PATH, sep=sep, nrows=0).columns)
    print("Columns:", cols)

    smiles_cols = [c for c in cols if "smiles" in c.lower()]
    if not smiles_cols:
        raise SystemExit("No SMILES column found. Send me: head -n 2 data/raw/coconut.csv | cut -c1-600")
    preferred = [c for c in smiles_cols if c.lower() == "canonical_smiles"]
    smiles_col = preferred[0] if preferred else smiles_cols[0]
    print("Using SMILES column:", smiles_col)

    keep = [smiles_col] + [c for c in cols if c != smiles_col
                           and any(h in c.lower() for h in KEEP_HINTS)]
    df = pd.read_csv(IN_PATH, sep=sep, usecols=keep, low_memory=False)
    print(f"Loaded {len(df):,} rows, keeping columns: {keep}")

    uniq = df[smiles_col].dropna().drop_duplicates()
    mapping = {s: standardize(s) for s in tqdm(uniq, desc="Standardizing")}
    df["smiles_std"] = df[smiles_col].map(mapping)
    df = df.dropna(subset=["smiles_std"]).drop_duplicates(subset=["smiles_std"])

    df.to_parquet(OUT_PATH, index=False)
    print(f"Saved {len(df):,} natural products")


if __name__ == "__main__":
    main()