import os
import sys
from pathlib import Path

import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw

ROOT = Path(__file__).resolve().parent.parent
os.chdir(ROOT)
sys.path.insert(0, str(ROOT / "src"))

from admet import get_admet_flags  # noqa: E402
from predict import predict_targets, _load  # noqa: E402
from utils import standardize  # noqa: E402

EXAMPLES = {
    "Caffeine": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "Capsaicin": "COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O",
    "Piperine": "O=C(/C=C/C=C/c1ccc2OCOc2c1)N1CCCCC1",
    "Resveratrol": "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1",
    "Quercetin": "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",
    "Curcumin": "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O",
}

st.set_page_config(page_title="NP Target Navigator", layout="wide")


@st.cache_resource(show_spinner="Loading ligand library (first run only)...")
def warm_up():
    _load()
    return True


warm_up()

st.title("NP Target Navigator")
st.caption("Similarity-based protein target hypotheses and drug-likeness flags for natural compounds.")

with st.sidebar:
    st.header("Settings")
    top_n = st.slider("Targets to show", 5, 20, 10)
    k = st.slider("Nearest known ligands used", 10, 100, 50, step=10)
    min_sim = st.slider("Minimum similarity", 0.2, 0.6, 0.3, step=0.05)
    hide_exact = st.checkbox(
        "Hide exact library matches (treat the compound as unseen)", value=True)
    st.markdown("---")
    st.caption("Knowledge base: human single-protein ChEMBL activities (pChEMBL >= 5). "
               "Natural product structures: COCONUT.")

choice = st.selectbox("Start from an example (optional)", ["Custom SMILES"] + list(EXAMPLES))
smiles = st.text_input("SMILES", value=EXAMPLES.get(choice, ""))

if st.button("Predict targets", type="primary"):
    if not smiles.strip():
        st.error("Enter a SMILES string.")
        st.stop()
    try:
        std = standardize(smiles.strip())
        if std is None:
            raise ValueError("Could not parse that SMILES string.")
        with st.spinner("Searching 900,000+ known ligands..."):
            res = predict_targets(std, top_n=top_n, k=k, min_sim=min_sim,
                                  exclude_exact=hide_exact)
        flags = get_admet_flags(std)
    except Exception as e:
        st.error(str(e))
        st.stop()

    max_sim = res.attrs.get("max_sim", 0.0)
    left, right = st.columns([1, 2])

    with left:
        st.subheader("Structure")
        st.image(Draw.MolToImage(Chem.MolFromSmiles(std), size=(380, 300)))
        st.subheader("Drug-likeness (rule-based)")
        c1, c2 = st.columns(2)
        c1.metric("Mol. weight", flags["mol_weight"])
        c2.metric("LogP", flags["logp"])
        c1.metric("QED", flags["qed"])
        c2.metric("TPSA", flags["tpsa"])
        c1.metric("Lipinski violations", flags["lipinski_violations"])
        c2.metric("Veber OK", "yes" if flags["veber_ok"] else "no")
        st.info(flags["summary"])
        if flags["pains_alerts"]:
            st.write("PAINS:", ", ".join(flags["pains_alerts"]))
        if flags["brenk_alerts"]:
            st.write("Structural alerts:", ", ".join(flags["brenk_alerts"]))

    with right:
        st.subheader("Predicted targets")
        if max_sim < 0.4:
            st.warning(f"Closest known ligand has similarity {max_sim:.2f}. "
                       "This compound is far from the training data, so treat these "
                       "predictions as low-confidence.")
        else:
            st.caption(f"Closest known ligand similarity: {max_sim:.2f}")
        if res.empty:
            st.write("No known ligands above the similarity threshold. Try lowering it.")
        else:
            st.dataframe(res.round(2), hide_index=True)
            chart = res.assign(label=res["target_name"].str.slice(0, 35)
                               + " (" + res["target_chembl_id"] + ")")
            st.bar_chart(chart.set_index("label")["confidence"])
            st.download_button("Download predictions (CSV)",
                               res.to_csv(index=False), "predicted_targets.csv")

    st.caption("Hypotheses for experimental follow-up, not confirmed activities. "
               "Confidence is the share of neighbour similarity supporting a target, "
               "not a probability. ADMET flags are rule-based.")