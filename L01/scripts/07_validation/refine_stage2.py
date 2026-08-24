# -*- coding: utf-8 -*-
# refine_stage2.py — multi-copy 선정 pose의 PoseBusters 결함 copy만 gnina 국소최소화로 정제.
#   5b_refine.py의 multi-copy 버전: posebusters.csv에서 valid=False인 (binder,copy)만 refine.
#   각 copy를 gnina --minimize(국소, iters 작게) → drift>max면 원본유지 → 그 copy heavy atom만 cif에 덮어씀.
#   나머지 copy·단백질·Zn 불변. refined cif per binder + stage2_mc_selection_refined.csv 생성.
#   boltz2 env(gemmi+rdkit) + singularity gnina.
import argparse, csv, os, subprocess
import gemmi
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds, rdMolAlign
from rdkit.Geometry import Point3D
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')

GNINA_SIF = os.environ.get("GNINA_SIF", "/path/to/casp17-ligand/models/gnina/gnina.sif")
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


def gnina_minimize(prot, lig_sdf, out_sdf, iters=10):
    cmd = ["singularity", "exec", "--bind", "/path/to", GNINA_SIF, "gnina",
           "--receptor", prot, "--ligand", lig_sdf, "--minimize",
           "--minimize_iters", str(iters), "--autobox_ligand", lig_sdf, "--out", out_sdf]
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r.stderr


def drift_and_coords(input_mol, refined_sdf):
    ref = next(Chem.SDMolSupplier(refined_sdf, removeHs=True))
    inp = Chem.RemoveHs(Chem.Mol(input_mol))
    try:
        d = rdMolAlign.CalcRMS(ref, inp)
    except Exception:
        d = 999.0
    conf = ref.GetConformer(); match = ref.GetSubstructMatch(inp)
    if len(match) == inp.GetNumAtoms():
        coords = [[conf.GetAtomPosition(j).x, conf.GetAtomPosition(j).y, conf.GetAtomPosition(j).z]
                  for j in match]
    else:
        coords = [[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
                  for i in range(ref.GetNumAtoms())]
    return d, coords


def organic_residues(model):
    out = []
    for ch in model:
        for r in ch:
            t = gemmi.find_tabulated_residue(r.name)
            if (t and (t.is_amino_acid() or t.is_nucleic_acid())) or r.name in WATER or r.name == "ZN":
                continue
            if any(a.element.name == "C" for a in r):
                out.append(r)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sel", required=True)
    ap.add_argument("--smiles", required=True)
    ap.add_argument("--posebusters", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-drift", type=float, default=2.0)
    ap.add_argument("--iters", type=int, default=10)
    ap.add_argument("--only", default="", help="쉼표구분 binder만(테스트용)")
    a = ap.parse_args()
    work = os.path.join(a.out, "work"); os.makedirs(work, exist_ok=True)

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
        reslist = organic_residues(st[0])
        pst = gemmi.read_structure(cif); pst.setup_entities()
        pst.remove_ligands_and_waters(); pst.remove_empty_chains()
        prot = os.path.join(work, f"{b}_prot.pdb"); pst.write_pdb(prot)
        changed = 0
        for i, res in enumerate(reslist):
            if i not in fail[b]:
                continue
            atoms = [(at.element.name, (at.pos.x, at.pos.y, at.pos.z))
                     for at in res if at.element.name != "H"]
            try:
                mol = build_mol(s, atoms)
                sdf = os.path.join(work, f"{b}_c{i}.sdf")
                w = Chem.SDWriter(sdf); w.write(mol); w.close()
                rsdf = os.path.join(work, f"{b}_c{i}_ref.sdf")
                err = gnina_minimize(prot, sdf, rsdf, a.iters)
                if not os.path.exists(rsdf):
                    rep.append([b, i, "", f"gnina_no_out:{err.strip()[-60:]}"]); continue
                d, coords = drift_and_coords(mol, rsdf)
                if d <= a.max_drift:
                    heavy = [at for at in res if at.element.name != "H"]
                    for k, at in enumerate(heavy):
                        if k < len(coords):
                            at.pos = gemmi.Position(*coords[k])
                    changed += 1; rep.append([b, i, round(d, 3), "refined"])
                else:
                    rep.append([b, i, round(d, 3), f"orig(drift>{a.max_drift})"])
            except Exception as e:
                rep.append([b, i, "", f"err:{str(e)[:50]}"])
        out_cif = os.path.join(a.out, f"{b}.cif")
        st.make_mmcif_document().write_file(out_cif)
        print(f"[{b}] {changed}/{len(fail[b])} copy refined -> {out_cif}")

    outsel = os.path.join(os.path.dirname(a.sel), "stage2_mc_selection_refined.csv")
    cols = list(selrows[0].keys())
    with open(outsel, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for row in selrows:
            b = row["binder"]; rc = os.path.join(a.out, f"{b}.cif")
            if b in fail and os.path.exists(rc):
                row["cif"] = rc
            w.writerow({k: row.get(k, "") for k in cols})
    with open(os.path.join(a.out, "refine_summary.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["binder", "copy", "drift", "used"]); w.writerows(rep)
    print(f"[refine] refined cif -> {a.out} , 새 선정 -> {outsel}")


if __name__ == "__main__":
    main()
