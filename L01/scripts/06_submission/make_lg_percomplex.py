# -*- coding: utf-8 -*-
# CASP17 LG stage-2 (per-complex, multi-copy + Zn) — Example 6.1 형식.
#   binder마다: 자기 수용체(PDB) + 유기 리간드 N copy(각 MDL) + Zn 이온(단일원자 MDL, M CHG)
#   -> ./L01/<id>LG<grp>_1
# 입력: stage2_mc_selection.csv(선정, cif 전체경로) + L01.smiles.stage2.csv(HTX명·SMILES·stoich)
# env: CASP_GROUP, CASP_AUTHOR, (선택) L01_BASE/L01_SEL/L01_CSV/L01_OUT
import gemmi, csv, os, sys, shlex
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds
from rdkit.Geometry import Point3D
from rdkit import RDLogger; RDLogger.DisableLog('rdApp.*')

BASE = os.environ.get("L01_BASE", "/path/to/casp17-ligand/users/USERNAME/targets/L01")
CON  = BASE + "/binders"
SEL  = os.environ.get("L01_SEL", CON + "/stage2_mc_selection.csv")
CSV  = os.environ.get("L01_CSV", CON + "/L01.smiles.stage2.csv")
OUTROOT = os.environ.get("L01_OUT", CON + "/stage2_submit/L01")
GROUP  = os.environ.get("CASP_GROUP", "YYY")
AUTHOR = os.environ.get("CASP_AUTHOR", "CHANGE-ME")
os.makedirs(OUTROOT, exist_ok=True)

# binder -> (HTX 코드, SMILES, stoichiometry)  from stage2 CSV
lig_info = {}
with open(CSV, newline="") as f:
    rd = csv.reader(f); next(rd, None)
    for row in rd:
        if not row or not row[0].strip().startswith("L01"):
            continue
        try: st = int(row[3])
        except (IndexError, ValueError): st = 1
        lig_info[row[0].strip()] = (row[1].strip(), row[2].strip(), st)


def norm_el(e):
    e = e.strip()
    return e[0].upper() + e[1:].lower() if e else e


def parse_atom_site(cif):
    rows = []; cols = []; ins = False
    for raw in open(cif):
        line = raw.strip()
        if not line:
            continue
        if line == "loop_":
            cols = []; ins = False; continue
        if raw.startswith("_atom_site."):
            cols.append(raw.strip().split(".")[-1]); ins = True; continue
        if not ins or not cols:
            continue
        if raw.startswith("_"):
            if rows: break
            continue
        p = shlex.split(raw, posix=False)
        if len(p) == len(cols):
            rows.append(dict(zip(cols, p)))
        elif rows:
            break
    return rows


def receptor_pdb_lines(rows):
    out = []; s = 1
    for r in rows:
        if r["group_PDB"] != "ATOM":
            continue
        el = norm_el(r["type_symbol"])
        if el == "H":
            continue
        nm = r["label_atom_id"].strip()[:4]
        nm = f" {nm:<3}" if len(el) == 1 and len(nm) < 4 else f"{nm:<4}"
        res = r["label_comp_id"][:3].upper()
        ch = (r.get("auth_asym_id") or r.get("label_asym_id") or "A")[:1]
        try: rs = int(float(r.get("auth_seq_id") or r.get("label_seq_id") or s))
        except Exception: rs = s
        x, y, z = float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"])
        b = float(r.get("B_iso_or_equiv", 50.0))
        ln = f"ATOM  {s:5d} {nm} {res:>3} {ch}{rs:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00{b:6.2f}          {el:>2}"
        out.append(ln[:80].ljust(80)); s += 1
    return out


def build_mol(smiles, atoms):
    """유기 copy 1개(heavy atoms)를 SMILES로 mol화(좌표 배정, 결합차수 복원)."""
    tmpl = Chem.RemoveHs(Chem.MolFromSmiles(smiles)); n = tmpl.GetNumAtoms()
    if len(atoms) != n:
        raise ValueError(f"atom {n} vs {len(atoms)}")
    same = all(norm_el(tmpl.GetAtomWithIdx(i).GetSymbol()) == norm_el(el)
               for i, (el, _) in enumerate(atoms))
    if same:
        m = Chem.Mol(tmpl); c = Chem.Conformer(n)
        for i, (_, xyz) in enumerate(atoms):
            c.SetAtomPosition(i, Point3D(*xyz))
        m.RemoveAllConformers(); m.AddConformer(c, assignId=True); return m
    rw = Chem.RWMol(); c = Chem.Conformer(n)
    for i, (el, xyz) in enumerate(atoms):
        rw.AddAtom(Chem.Atom(norm_el(el))); c.SetAtomPosition(i, Point3D(*xyz))
    m = rw.GetMol(); m.AddConformer(c, assignId=True)
    rdDetermineBonds.DetermineConnectivity(m)
    return AllChem.AssignBondOrdersFromTemplate(tmpl, m)


def mol_lines(m):
    L = [l.rstrip() for l in Chem.MolToMolBlock(m).splitlines()]
    while L and not L[-1]:
        L.pop()
    return L


def ion_mol_lines(el, xyz, charge=2):
    """단일 금속이온 MDL (M CHG 포함). Zn=+2."""
    rw = Chem.RWMol(); a = Chem.Atom(norm_el(el)); a.SetFormalCharge(charge)
    rw.AddAtom(a); m = rw.GetMol()
    c = Chem.Conformer(1); c.SetAtomPosition(0, Point3D(*xyz)); m.AddConformer(c, assignId=True)
    return mol_lines(m)


def split_ligands(st):
    """비단백질 잔기를 organic copy들 + 이온(무탄소)으로 분리.
    반환: (copies=[[(el,xyz)...], ...], ions=[(el,xyz), ...])."""
    copies = []; ions = []
    for ch in st[0]:
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            if info and info.is_amino_acid():
                continue
            if res.name in ("HOH", "WAT", "DOD", "H2O"):
                continue
            atoms = [(at.element.name, (at.pos.x, at.pos.y, at.pos.z))
                     for at in res if at.element.name != "H"]
            if not atoms:
                continue
            if any(el == "C" for el, _ in atoms):
                copies.append(atoms)                     # 유기 copy
            else:
                ions += atoms                            # 금속/이온(Zn 등)
    return copies, ions


rows = list(csv.DictReader(open(SEL, encoding="utf-8"))); ok = 0; fail = []; warn = []
for r in rows:
    b = r["binder"]; cif = r.get("cif", "")
    if not cif or not os.path.exists(cif):
        fail.append((b, "cif")); continue
    htx, smi, stoich = lig_info.get(b, ("", "", 1))
    if not smi:
        fail.append((b, "smiles")); continue
    try:
        st = gemmi.read_structure(cif); st.setup_entities()
    except Exception as e:
        fail.append((b, f"read:{e}")); continue
    copies, ions = split_ligands(st)
    if not copies:
        fail.append((b, "no-organic")); continue
    if len(copies) != stoich:
        warn.append((b, f"copies {len(copies)}!=stoich {stoich}"))
    rec = receptor_pdb_lines(parse_atom_site(cif))
    try:
        ls = f"{max(0.0, min(1.0, float(r.get('ligand_iptm', '0.5')))):.3f}"
    except Exception:
        ls = "0.500"

    body = ["PFRMAT LG", f"TARGET {b}", f"AUTHOR {AUTHOR}",
            "METHOD Consensus pose from co-folding (Boltz2, Protenix, AlphaFold3), pocket-based selection.",
            "MODEL 1"]
    body += rec + ["TER"]
    n = 0; ok_all = True
    for atoms in copies:                                 # 유기 copy들: LIGAND 1..N
        try:
            m = build_mol(smi, atoms)
        except Exception as e:
            fail.append((b, f"copy:{e}")); ok_all = False; break
        n += 1
        body += [f"LIGAND {n} {htx}", f"LSCORE {ls}"] + mol_lines(m)
    if not ok_all:
        continue
    for el, xyz in ions:                                 # 이온(Zn): LIGAND N+1..
        n += 1
        body += [f"LIGAND {n} {norm_el(el).upper()}", f"LSCORE {ls}"] + ion_mol_lines(el, xyz)
    body += ["END"]

    open(os.path.join(OUTROOT, f"{b}LG{GROUP}_1"), "w").write("\n".join(body) + "\n")
    ok += 1

print(f"생성 {ok}개 -> {OUTROOT}/<id>LG{GROUP}_1")
if warn:
    print(f"[주의] copy수≠stoich {len(warn)}개: {warn[:10]}")
if fail:
    print(f"[실패] {len(fail)}개: {fail[:10]}")

# ── 실행 ── (rdkit+gemmi 있는 env; 예: ligand_select 또는 boltz2)
#   CASP_GROUP=YYY CASP_AUTHOR=CHANGE-ME python make_lg_percomplex.py
#   그다음:  cd $L01_OUT/..; tar -czf L01LG_YYY.tgz ./L01
