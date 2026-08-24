# -*- coding: utf-8 -*-
# refine_mmff_stage2.py — gnina refine 후에도 남은 covalent 기하 결함 copy를 RDKit MMFF로 교정.
#   ring 평탄·bond length 실패 대상. heavy atom을 원위치에 제약(0.5Å)한 채 MMFF 최소화 →
#   결합길이·각도·고리 pucker 이상화(위치는 거의 유지). protein clash는 못 고침(→ A로).
#   입력 sel = gnina-refined 선정(stage2_mc_selection_refined.csv), posebusters = 그 재검 결과.
#   boltz2 env(gemmi+rdkit)만. GPU/gnina 불필요.
import argparse, csv, os, math
import gemmi
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds
from rdkit.Geometry import Point3D
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')
WATER = ("HOH", "WAT", "DOD", "H2O")


def build_mol(smiles, atoms):
    tpl = Chem.RemoveHs(Chem.MolFromSmiles(smiles)); n = tpl.GetNumAtoms()
    if len(atoms) != n:
        raise ValueError(f"atom {n} vs {len(atoms)}")
    same = all(tpl.GetAtomWithIdx(i).GetSymbol().capitalize() == el.capitalize()
               for i, (el, _) in enumerate(atoms))
    if same:
        m = Chem.Mol(tpl); c = Chem.Conformer(n)
        for i, (_, xyz) in enumerate(atoms):
            c.SetAtomPosition(i, Point3D(*[float(v) for v in xyz]))
        m.RemoveAllConformers(); m.AddConformer(c, assignId=True); return m
    rw = Chem.RWMol(); c = Chem.Conformer(n)
    for i, (el, xyz) in enumerate(atoms):
        rw.AddAtom(Chem.Atom(el.capitalize())); c.SetAtomPosition(i, Point3D(*[float(v) for v in xyz]))
    m = rw.GetMol(); m.AddConformer(c, assignId=True)
    rdDetermineBonds.DetermineConnectivity(m)
    return AllChem.AssignBondOrdersFromTemplate(tpl, m)


def mmff_refine(mol, maxdispl=0.5, fc=100.0, iters=1000):
    """heavy atom을 원위치±maxdispl로 제약한 MMFF 최소화. 실패 시 None."""
    mh = Chem.AddHs(mol, addCoords=True)
    mp = AllChem.MMFFGetMoleculeProperties(mh)
    if mp is None:
        return None
    ff = AllChem.MMFFGetMoleculeForceField(mh, mp)
    if ff is None:
        return None
    for at in mh.GetAtoms():
        if at.GetAtomicNum() > 1:
            ff.MMFFAddPositionConstraint(at.GetIdx(), maxdispl, fc)
    ff.Minimize(maxIts=iters)
    return Chem.RemoveHs(mh)   # heavy atom 순서 = 입력 순서 유지


def organics(m):
    o = []
    for ch in m:
        for r in ch:
            t = gemmi.find_tabulated_residue(r.name)
            if (t and (t.is_amino_acid() or t.is_nucleic_acid())) or r.name in WATER or r.name == "ZN":
                continue
            if any(a.element.name == "C" for a in r):
                o.append(r)
    return o


def main():
    ap = argparse.ArgumentParser()
    for x in ("--sel", "--smiles", "--posebusters", "--out"):
        ap.add_argument(x, required=True)
    ap.add_argument("--max-drift", type=float, default=1.0)
    ap.add_argument("--only", default="")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    smi = {}
    with open(a.smiles, newline="") as f:
        rd = csv.reader(f); next(rd, None)
        for row in rd:
            if row and row[0].startswith("L01"):
                smi[row[0]] = row[2]

    fail = {}
    for r in csv.DictReader(open(a.posebusters)):
        if str(r.get("valid")) == "False":
            fail.setdefault(r["pocket_id"], set()).add(int(r["ligand_cluster_id"]))
    only = set(x.strip() for x in a.only.split(",") if x.strip())

    selrows = list(csv.DictReader(open(a.sel, encoding="utf-8")))
    rep = []
    for row in selrows:
        b = row["binder"]; cif = row.get("cif", "")
        if b not in fail or not cif or not os.path.exists(cif):
            continue
        if only and b not in only:
            continue
        s = smi.get(b, "")
        st = gemmi.read_structure(cif); st.setup_entities()
        reslist = organics(st[0])
        changed = 0
        for i, res in enumerate(reslist):
            if i not in fail[b]:
                continue
            atoms = [(at.element.name, (at.pos.x, at.pos.y, at.pos.z))
                     for at in res if at.element.name != "H"]
            try:
                mol = build_mol(s, atoms)
                rmol = mmff_refine(mol)
                if rmol is None:
                    rep.append([b, i, "", "mmff_unparam"]); continue
                conf = rmol.GetConformer()
                coords = [[conf.GetAtomPosition(k).x, conf.GetAtomPosition(k).y, conf.GetAtomPosition(k).z]
                          for k in range(rmol.GetNumAtoms())]
                d = math.sqrt(sum((coords[k][0]-atoms[k][1][0])**2 + (coords[k][1]-atoms[k][1][1])**2
                                  + (coords[k][2]-atoms[k][1][2])**2 for k in range(len(atoms))) / len(atoms))
                if d <= a.max_drift:
                    heavy = [at for at in res if at.element.name != "H"]
                    for k, at in enumerate(heavy):
                        if k < len(coords):
                            at.pos = gemmi.Position(*coords[k])
                    changed += 1; rep.append([b, i, round(d, 3), "mmff"])
                else:
                    rep.append([b, i, round(d, 3), f"orig(drift>{a.max_drift})"])
            except Exception as e:
                rep.append([b, i, "", f"err:{str(e)[:50]}"])
        st.make_mmcif_document().write_file(os.path.join(a.out, f"{b}.cif"))
        print(f"[{b}] {changed}/{len(fail[b])} copy MMFF 교정")

    outsel = os.path.join(os.path.dirname(a.sel), "stage2_mc_selection_mmff.csv")
    cols = list(selrows[0].keys())
    with open(outsel, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for row in selrows:
            b = row["binder"]; rc = os.path.join(a.out, f"{b}.cif")
            if b in fail and os.path.exists(rc):
                row["cif"] = rc
            w.writerow({k: row.get(k, "") for k in cols})
    with open(os.path.join(a.out, "refine_mmff_summary.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["binder", "copy", "drift", "used"]); w.writerows(rep)
    print(f"[mmff] -> {a.out} , 새 선정 -> {outsel}")


if __name__ == "__main__":
    main()
