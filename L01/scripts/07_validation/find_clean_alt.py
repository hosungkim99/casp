# -*- coding: utf-8 -*-
# find_clean_alt.py — 지정 binder에 대해 후보 구조 중 (posebusters 전 copy valid + exosite<near)인
#   대체 pose를 순위대로 탐색, 첫 통과 후보 채택. gnina/MMFF로 못 고친 clash binder의 A(교체)용.
#   후보 = stage2_mc_candidates.csv(cif_full 컬럼 필요, --top 크게 재생성).
#   boltz2 env(gemmi+rdkit) + 4c posebusters(casp_eval, subprocess 호출).
import argparse, csv, os, subprocess, math
import gemmi
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds
from rdkit.Geometry import Point3D
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')

DRUG = {"7X9", "7WK", "6JP"}
REFIDS = ["7POL", "7POO", "7POQ", "7POU"]


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


def extract(cif):
    st = gemmi.read_structure(cif); st.setup_entities(); copies = []
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


def load_ref(refdir, rid):
    for ext in (".pdb", ".cif"):
        p = os.path.join(refdir, rid + ext)
        if os.path.exists(p):
            st = gemmi.read_structure(p); st.setup_entities()
            drug = [at.pos for ch in st[0] for res in ch if res.name in DRUG
                    for at in res if at.element.name != "H"]
            if drug:
                return st[0][0].get_polymer(), drug
    return None


def exosite_min(cif, refs):
    st = gemmi.read_structure(cif); st.setup_entities(); op = st[0][0].get_polymer()
    lig = []
    for ch in st[0]:
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            if (info and info.is_amino_acid()) or res.name in ("HOH", "WAT", "DOD", "H2O", "ZN"):
                continue
            if any(at.element.name == "C" for at in res):
                lig += [at.pos for at in res if at.element.name != "H"]
    if not lig:
        return 999.0
    best = 1e9
    for rp, drug in refs:
        sup = gemmi.calculate_superposition(op, rp, gemmi.PolymerType.PeptideL, gemmi.SupSelect.CaP)
        dt = [sup.transform.apply(dp) for dp in drug]
        for lp in lig:
            for v in dt:
                d = math.sqrt((lp.x-v.x)**2 + (lp.y-v.y)**2 + (lp.z-v.z)**2)
                if d < best:
                    best = d
    return best


def prep(cif, smiles, wd, tag):
    copies, prot = extract(cif); os.makedirs(wd, exist_ok=True)
    pdb = os.path.join(wd, f"{tag}_prot.pdb"); prot.write_pdb(pdb); rows = []
    for i, atoms in enumerate(copies):
        sdf = os.path.join(wd, f"{tag}_c{i}.sdf")
        try:
            m = build_mol(smiles, atoms); w = Chem.SDWriter(sdf); w.write(m); w.close()
        except Exception:
            sdf = ""
        rows.append({"pocket_id": tag, "ligand_cluster_id": i, "member_rank": 1,
                     "rep_cif": cif, "sdf": sdf, "protein_pdb": pdb, "note": ""})
    return rows


def run_pb(rows, wd, pbpy, pbsc):
    man = os.path.join(wd, "manifest.csv")
    with open(man, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["pocket_id", "ligand_cluster_id", "member_rank",
                                          "rep_cif", "sdf", "protein_pdb", "note"])
        w.writeheader(); w.writerows(rows)
    out = os.path.join(wd, "pb.csv")
    subprocess.run([pbpy, pbsc, "--manifest", man, "--out", out], capture_output=True, text=True)
    if not os.path.exists(out):
        return None
    return list(csv.DictReader(open(out)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--smiles", required=True)
    ap.add_argument("--refdir", required=True)
    ap.add_argument("--binders", required=True, help="쉼표구분")
    ap.add_argument("--out", required=True)
    ap.add_argument("--pb-python", required=True)
    ap.add_argument("--pb-script", required=True)
    ap.add_argument("--top", type=int, default=30)
    ap.add_argument("--near", type=float, default=8.0)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    refs = [r for r in (load_ref(a.refdir, x) for x in REFIDS) if r]
    if not refs:
        raise SystemExit("[에러] ref 약물구조 없음")
    smi = {}
    with open(a.smiles, newline="") as f:
        rd = csv.reader(f); next(rd, None)
        for row in rd:
            if row and row[0].startswith("L01"):
                smi[row[0]] = row[2]
    cands = {}
    for r in csv.DictReader(open(a.candidates)):
        cands.setdefault(r["binder"], []).append(r)
    for b in cands:
        cands[b].sort(key=lambda r: int(r["submit_rank"]))

    targets = [x.strip() for x in a.binders.split(",") if x.strip()]
    picks = []
    for b in targets:
        chosen = None
        for c in cands.get(b, [])[:a.top]:
            cif = c.get("cif_full", "")
            rk = c["submit_rank"]
            if not cif or not os.path.exists(cif):
                continue
            wd = os.path.join(a.out, f"{b}_r{rk}")
            rows = prep(cif, smi.get(b, ""), wd, f"{b}_r{rk}")
            pb = run_pb(rows, wd, a.pb_python, a.pb_script)
            allvalid = bool(pb) and len(pb) > 0 and all(str(x["valid"]) == "True" for x in pb)
            exo = exosite_min(cif, refs) if allvalid else None
            print(f"  {b} r{rk}: valid={allvalid} exosite={round(exo,2) if exo is not None else '-'}")
            if allvalid and exo is not None and exo <= a.near:
                chosen = (rk, cif, round(exo, 2)); break
        if chosen:
            picks.append({"binder": b, "rank": chosen[0], "cif": chosen[1],
                          "exosite": chosen[2], "status": "found"})
        else:
            picks.append({"binder": b, "rank": "", "cif": "", "exosite": "", "status": "none"})
        print(f"[{b}] -> {picks[-1]['status']} (rank {picks[-1]['rank']})")

    with open(os.path.join(a.out, "clean_alt.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["binder", "rank", "cif", "exosite", "status"])
        w.writeheader(); w.writerows(picks)
    print(f"[find_clean_alt] -> {os.path.join(a.out, 'clean_alt.csv')}")


if __name__ == "__main__":
    main()
