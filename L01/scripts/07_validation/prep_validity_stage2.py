# -*- coding: utf-8 -*-
# prep_validity_stage2.py — multi-copy 선정 pose → PoseBusters 입력(copy별 SDF + 단백질 PDB) + manifest.
#   boltz2 env(gemmi+rdkit). 이후 boost/4c_posebusters.py(casp_eval env)가 manifest로 검사.
#   4b_prep_validity.py의 multi-copy·stage2 버전: ligand_clusters 대신 stage2_mc_selection.csv,
#   copy 1개가 아니라 유기 copy 전부를 각각 SDF로. manifest 컬럼은 4c가 읽는 형식 그대로.
import argparse, csv, os
import gemmi
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds
from rdkit.Geometry import Point3D
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')


def build_mol(smiles, atoms):
    tmpl = Chem.RemoveHs(Chem.MolFromSmiles(smiles)); n = tmpl.GetNumAtoms()
    if len(atoms) != n:
        raise ValueError(f"atom {n} vs {len(atoms)}")
    same = all(tmpl.GetAtomWithIdx(i).GetSymbol().capitalize() == el.capitalize()
               for i, (el, _) in enumerate(atoms))
    if same:
        m = Chem.Mol(tmpl); c = Chem.Conformer(n)
        for i, (_, xyz) in enumerate(atoms):
            c.SetAtomPosition(i, Point3D(*[float(v) for v in xyz]))
        m.RemoveAllConformers(); m.AddConformer(c, assignId=True); return m
    rw = Chem.RWMol(); c = Chem.Conformer(n)
    for i, (el, xyz) in enumerate(atoms):
        rw.AddAtom(Chem.Atom(el.capitalize())); c.SetAtomPosition(i, Point3D(*[float(v) for v in xyz]))
    m = rw.GetMol(); m.AddConformer(c, assignId=True)
    rdDetermineBonds.DetermineConnectivity(m)
    return AllChem.AssignBondOrdersFromTemplate(tmpl, m)


def extract(cif):
    """유기 copy들[(el,(x,y,z))...] + 단백질 Structure(리간드·물·이온 제거)."""
    st = gemmi.read_structure(cif); st.setup_entities()
    copies = []
    for ch in st[0]:
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            if (info and info.is_amino_acid()) or res.name in ("HOH", "WAT", "DOD", "H2O", "ZN"):
                continue
            atoms = [(at.element.name, (at.pos.x, at.pos.y, at.pos.z))
                     for at in res if at.element.name != "H"]
            if any(e == "C" for e, _ in atoms):
                copies.append(atoms)
    prot = gemmi.read_structure(cif); prot.setup_entities()
    prot.remove_ligands_and_waters(); prot.remove_empty_chains()
    return copies, prot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sel", required=True, help="stage2_mc_selection.csv")
    ap.add_argument("--smiles", required=True, help="L01.smiles.stage2.csv")
    ap.add_argument("--out", required=True, help="posebusters 폴더")
    a = ap.parse_args()
    inp = os.path.join(a.out, "inputs"); os.makedirs(inp, exist_ok=True)

    smi = {}
    with open(a.smiles, newline="") as f:
        rd = csv.reader(f); next(rd, None)
        for row in rd:
            if row and row[0].startswith("L01"):
                smi[row[0]] = row[2]

    man = []
    for r in csv.DictReader(open(a.sel, encoding="utf-8")):
        b = r["binder"]; cif = r.get("cif", "")
        if not cif or not os.path.exists(cif):
            continue
        s = smi.get(b, "")
        copies, prot = extract(cif)
        pdb = os.path.join(inp, f"{b}_protein.pdb"); prot.write_pdb(pdb)
        for i, atoms in enumerate(copies):
            tag = f"{b}_c{i}"; sdf = ""; note = ""
            if not s:
                note = "no_smiles"
            else:
                try:
                    m = build_mol(s, atoms)
                    sdf = os.path.join(inp, tag + ".sdf")
                    w = Chem.SDWriter(sdf); w.write(m); w.close()
                except Exception as e:
                    sdf = ""; note = f"prep_err:{str(e)[:50]}"
            man.append({"pocket_id": b, "ligand_cluster_id": i, "member_rank": 1,
                        "rep_cif": cif, "sdf": sdf, "protein_pdb": pdb, "note": note})
            print(f"  {tag}: sdf={'O' if sdf else 'X'} {note}")

    mcsv = os.path.join(a.out, "manifest.csv")
    with open(mcsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["pocket_id", "ligand_cluster_id", "member_rank",
                                          "rep_cif", "sdf", "protein_pdb", "note"])
        w.writeheader(); w.writerows(man)
    print(f"[prep] {len(man)}개 copy 준비 (binder별 protein.pdb + copy별 sdf) -> {mcsv}")


if __name__ == "__main__":
    main()
