from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, QED, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

from utils import standardize


def _make_catalog(which):
    params = FilterCatalogParams()
    params.AddCatalog(which)
    return FilterCatalog(params)


_PAINS = _make_catalog(FilterCatalogParams.FilterCatalogs.PAINS)
_BRENK = _make_catalog(FilterCatalogParams.FilterCatalogs.BRENK)


def _alerts(catalog, mol):
    return sorted({m.GetDescription() for m in catalog.GetMatches(mol)})


def get_admet_flags(smiles):
    """Rule-based drug-likeness and structural-alert summary for one compound."""
    std = standardize(smiles)
    if std is None:
        raise ValueError("Invalid SMILES")
    mol = Chem.MolFromSmiles(std)

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)

    lipinski_violations = int(mw > 500) + int(logp > 5) + int(hbd > 5) + int(hba > 10)
    veber_ok = rotb <= 10 and tpsa <= 140
    pains = _alerts(_PAINS, mol)
    brenk = _alerts(_BRENK, mol)

    notes = []
    if lipinski_violations > 1:
        notes.append(f"{lipinski_violations} Lipinski violations (poor oral absorption likely)")
    if not veber_ok:
        notes.append("fails Veber rules (TPSA/rotatable bonds)")
    if pains:
        notes.append("PAINS alert (may give false positives in assays)")
    if brenk:
        notes.append("structural alerts present (see list below; not necessarily toxic)")

    return {
        "smiles_std": std,
        "mol_weight": round(mw, 1),
        "logp": round(logp, 2),
        "hbd": hbd,
        "hba": hba,
        "tpsa": round(tpsa, 1),
        "rotatable_bonds": rotb,
        "fraction_csp3": round(rdMolDescriptors.CalcFractionCSP3(mol), 2),
        "stereocenters": len(Chem.FindMolChiralCenters(mol, includeUnassigned=True)),
        "qed": round(QED.qed(mol), 2),
        "lipinski_violations": lipinski_violations,
        "veber_ok": veber_ok,
        "pains_alerts": pains,
        "brenk_alerts": brenk,
        "summary": "; ".join(notes) if notes else "no major drug-likeness flags",
    }


if __name__ == "__main__":
    import json
    import sys
    smi = sys.argv[1] if len(sys.argv) > 1 else "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O"
    print(json.dumps(get_admet_flags(smi), indent=2))