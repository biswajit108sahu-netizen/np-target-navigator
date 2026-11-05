from admet import get_admet_flags
from predict import predict_targets

COMPOUNDS = {
    "Caffeine": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "Capsaicin": "COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O",
    "Piperine": "O=C(/C=C/C=C/c1ccc2OCOc2c1)N1CCCCC1",
    "Resveratrol": "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1",
    "Quercetin": "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12",
    "Curcumin": "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O",
}


def main():
    lines = []
    for name, smi in COMPOUNDS.items():
        res = predict_targets(smi, top_n=8, exclude_exact=True)
        flags = get_admet_flags(smi)
        lines.append(f"=== {name} ===")
        lines.append(f"Closest non-identical library compound: {res.attrs['max_sim']:.2f}")
        lines.append(f"Drug-likeness: {flags['summary']} (QED {flags['qed']})")
        lines.append(res[["target_name", "score", "confidence", "n_neighbors"]]
                     .round(2).to_string(index=False))
        lines.append("")
    text = "\n".join(lines)
    print(text)
    with open("docs/case_studies_output.txt", "w") as f:
        f.write(text)


if __name__ == "__main__":
    main()