import sqlite3

import pandas as pd
from tqdm import tqdm

from utils import standardize

DB_PATH = "data/raw/chembl.db"
OUT_PATH = "data/processed/activities_clean.parquet"
PCHEMBL_MIN = 5.0  # about IC50/Ki <= 10 uM

QUERY = f"""
SELECT cs.canonical_smiles,
       td.chembl_id  AS target_chembl_id,
       td.pref_name  AS target_name,
       act.standard_type,
       act.pchembl_value
FROM activities act
JOIN assays a               ON act.assay_id = a.assay_id
JOIN target_dictionary td   ON a.tid = td.tid
JOIN compound_structures cs ON act.molregno = cs.molregno
WHERE act.pchembl_value >= {PCHEMBL_MIN}
  AND act.standard_relation = '='
  AND td.target_type = 'SINGLE PROTEIN'
  AND td.organism = 'Homo sapiens'
  AND a.confidence_score >= 8
"""


def main():
    con = sqlite3.connect(DB_PATH)
    print("Running SQL query (can take a few minutes)...")
    df = pd.read_sql_query(QUERY, con)
    con.close()
    print(f"Rows after SQL filter: {len(df):,}")

    uniq = df["canonical_smiles"].drop_duplicates()
    mapping = {s: standardize(s) for s in tqdm(uniq, desc="Standardizing")}
    df["smiles_std"] = df["canonical_smiles"].map(mapping)
    df = df.dropna(subset=["smiles_std"])
    df = df.drop(columns=["canonical_smiles"])
    df = df.drop_duplicates(subset=["smiles_std", "target_chembl_id"])

    df.to_parquet(OUT_PATH, index=False)
    print(f"Saved {len(df):,} compound-target pairs")
    print(f"Unique compounds: {df['smiles_std'].nunique():,}")
    print(f"Unique targets:   {df['target_chembl_id'].nunique():,}")


if __name__ == "__main__":
    main()