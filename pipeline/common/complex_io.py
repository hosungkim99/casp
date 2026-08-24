#!/usr/bin/env python3
"""
complex_io.py - 범용 단백질-리간드 복합체 파서 (gemmi).

타겟 종류·체인 수·리간드 개수에 무관하게 동작하도록, 모든 스텝이 공유하는 공용 파서.
- 단백질: 모든 polymer(아미노산) 체인의 Cα
- 리간드: 아미노산/핵산/물이 아닌 모든 잔기를 'copy' 단위로 검출 (개수 무관)
- 리간드 좌표를 canonical 순서로 concat → 구조 간 RMSD/클러스터링 가능
- 리간드 ↔ SMILES를 원소 조성으로 매칭 (ligand.tsv에 SMILES가 여러 개여도)

주의(범용성의 한계): 동일 조성 리간드가 '여러 copy'일 때 copy 간 대응은
seqid 순서로 고정한다. 진짜 교환대칭(symmetry) 최적 배정은 하지 않음(보수적).
"""
from collections import Counter
import gemmi
import os 

WATER = {"HOH", "WAT", "DOD", "H2O"}

def reference_cif(rows, ref_file=None):
    """
    정렬 기준 구조의 cif 경로를 고른다 (모든 스텝이 같은 기준을 공유하도록).
    - ref_file(= 1_protein_cluster가 쓴 'PC공간 중심 구조')이 있으면 그것을 사용.
    - 없으면 rows[0](rank 1위)로 폴백 (하위호환/단독실행).
    중심 구조 기준은 정렬 왜곡이 적고, 모델간 비교 불가한 ranking_score 편향을 피한다.
    """
    if ref_file and os.path.exists(ref_file):
        p = open(ref_file).read().strip()
        if p and os.path.exists(p):
            return p
    return rows[0]["cif"]


def parse_complex(path):
    """
    반환:
      prot_ca : dict[(chain, resnum)] -> [x,y,z]   (모든 단백질 체인의 Cα)
      ligands : list of dict {key:(chain,resname,seqid), resname, atoms:[(element,[x,y,z]),...]}
    """
    st = gemmi.read_structure(path)
    if len(st) == 0:
        return {}, []
    m = st[0]
    prot_ca, ligands = {}, []
    for ch in m:
        for res in ch:
            t = gemmi.find_tabulated_residue(res.name)
            is_aa = bool(t) and t.is_amino_acid()
            is_na = bool(t) and t.is_nucleic_acid()
            if is_aa:
                for a in res:
                    if a.name == "CA":
                        prot_ca[(ch.name, res.seqid.num)] = [a.pos.x, a.pos.y, a.pos.z]
                        break
            elif is_na or res.name in WATER:
                continue                       # 핵산·물은 리간드 아님
            else:
                atoms = [(a.element.name, [a.pos.x, a.pos.y, a.pos.z])
                         for a in res if a.element.name != "H"]
                if atoms:
                    ligands.append({"key": (ch.name, res.name, res.seqid.num),
                                    "resname": res.name, "atoms": atoms})
    return prot_ca, ligands


def _canon_order(ligands):
    """리간드 copy를 (resname, chain, seqid)로 정렬한 인덱스 → 구조 간 일관 순서."""
    return sorted(range(len(ligands)),
                  key=lambda i: (ligands[i]["resname"],
                                 ligands[i]["key"][0], ligands[i]["key"][2]))


# ── organic 리간드 vs 금속/이온 cofactor 분리 (multi-copy + Zn 대응) ──────────
# CASP L01/L02 stage2: 구조마다 유기 리간드 N copy + 금속이온(Zn/SFG 등)이 들어있다.
# 금속이온을 리간드로 취급하면 centroid/서명/도킹이 오염되므로 geometry에서 분리한다.
# 판정: heavy atom에 탄소가 하나도 없으면 이온/금속 cofactor로 본다(3_의 템플릿 이온제외와 동일 관례).

def is_organic(lg):
    """리간드 copy가 유기물인지(탄소 1개 이상). 단원자 금속이온(ZN 등)=False."""
    return any(el == "C" for el, _ in lg["atoms"])


def organic_ligands(ligands):
    """유기 리간드 copy만(금속/이온 cofactor 제외). geometry·도킹은 이걸 써야 정확."""
    return [lg for lg in ligands if is_organic(lg)]


def ion_ligands(ligands):
    """금속/이온 cofactor copy만(탄소 없는 잔기, 예: ZN). 출력엔 남기되 geometry엔 미포함."""
    return [lg for lg in ligands if not is_organic(lg)]


def copy_centroids(ligands, organic_only=True):
    """리간드 copy '각각'의 무게중심 리스트 [[x,y,z],...].
    multi-copy를 하나로 뭉개지 않고 copy별 점을 반환(포켓 후보 = copy별 자리).
    organic_only면 금속/이온 제외."""
    src = organic_ligands(ligands) if organic_only else ligands
    out = []
    for lg in src:
        pts = [xyz for _, xyz in lg["atoms"]]
        n = len(pts)
        out.append([sum(p[k] for p in pts) / n for k in range(3)])
    return out


def ligand_concat(ligands):
    """
    유기 리간드 heavy atom을 canonical 순서로 이어붙임(금속/이온 제외).
    반환: (elements:list[str], coords:list[[x,y,z]])
    elements는 구조 간 '서명' — 같은 타겟이면 모든 구조에서 동일해야 함.
    """
    ligs = organic_ligands(ligands)          # 금속/이온은 서명·좌표에서 제외
    elems, coords = [], []
    for i in _canon_order(ligs):
        for el, xyz in ligs[i]["atoms"]:
            elems.append(el); coords.append(xyz)
    return elems, coords


def ligand_centroid(ligands):
    """유기 리간드 전체 heavy atom의 무게중심 [x,y,z] (금속/이온 제외; 없으면 None).
    주의: multi-copy면 모든 copy를 한 점으로 평균하므로, 포켓 후보엔 copy_centroids()를 써라."""
    pts = [xyz for lg in organic_ligands(ligands) for _, xyz in lg["atoms"]]
    if not pts:
        return None
    n = len(pts)
    return [sum(p[k] for p in pts) / n for k in range(3)]


def composition(atoms_or_elems):
    """원소 multiset(Counter). 리간드↔SMILES 조성 매칭용."""
    if atoms_or_elems and isinstance(atoms_or_elems[0], tuple):
        elems = [el for el, _ in atoms_or_elems]
    else:
        elems = list(atoms_or_elems)
    return Counter(e.capitalize() for e in elems)


def smiles_composition(smiles):
    """SMILES의 heavy-atom 원소 조성(Counter). rdkit 필요."""
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.RemoveHs(mol)
    return Counter(a.GetSymbol() for a in mol.GetAtoms())


def match_ligand_to_smiles(ligand, smiles_list):
    """
    리간드 한 copy를 후보 SMILES들과 원소 조성으로 매칭.
    smiles_list: [(id, name, smiles), ...]
    반환: (id, name, smiles) 또는 None.
    - 후보 1개면 그대로 사용(단일 리간드 타겟).
    - 여러 개면 조성이 정확히 일치하는 것을 선택.
    """
    if len(smiles_list) == 1:
        return smiles_list[0]
    lig_comp = composition(ligand["atoms"])
    for sid, name, smi in smiles_list:
        if smiles_composition(smi) == lig_comp:
            return (sid, name, smi)
    return None


def read_ligand_tsv(path):
    """ligand.tsv → [(ID, Name, SMILES), ...] (행 여러 개 = 리간드 여러 종)."""
    import csv
    out = []
    with open(path, newline="") as h:
        for r in csv.DictReader(h, delimiter="\t"):
            if r.get("SMILES"):
                out.append((r.get("ID", ""), r.get("Name", ""), r["SMILES"]))
    return out


def write_ligands_pdb(ligands, pdb_out):
    """유기 리간드 copy를 한 PDB(HETATM, 체인 L)로 저장. gnina/도킹용.
    금속/이온(ZN 등)은 도킹 리간드에서 제외 — 리간드로 도킹하면 affinity가 왜곡됨."""
    lines = []
    serial = 1
    for ridx, lg in enumerate(organic_ligands(ligands), 1):
        rn = lg["resname"][:3].upper()
        for el, (x, y, z) in lg["atoms"]:
            nm = f"{el:<2}"
            lines.append(
                f"HETATM{serial:5d} {nm:<4} {rn:>3} L{ridx:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2}")
            serial += 1
    lines.append("END")
    open(pdb_out, "w").write("\n".join(lines) + "\n")
    return serial - 1

def write_receptor_pdb(cif, pdb_out):
    """cif에서 리간드·물 제거한 단백질만 PDB로 저장 (gnina 수용체/도킹 입력용).
    알고리즘: setup_entities → remove_ligands_and_waters → remove_empty_chains → write_pdb.
    (3·4b·5·5b·stage1 스텝이 공통으로 쓰던 코드를 한 곳으로 통일)"""
    st = gemmi.read_structure(cif); st.setup_entities()
    st.remove_ligands_and_waters(); st.remove_empty_chains(); st.write_pdb(pdb_out)


# ── 기하 캐시 (중복 파싱 제거) ──────────────────────────────────────────────
# 1·2·6·8 스텝은 같은 수천 개 cif를 각자 gemmi로 재파싱한다(gemmi 파싱이 병목).
# 0b_cache_geometry.py 가 한 번만 파싱해 만든 캐시를 이 함수들로 공유 → 재파싱 제거.
# 캐시 한 항목: {pk:[(chain,resnum)...], px:(M,3), le:[elem...], lx:(N,3)} (좌표만; 원본보다 작음)

def load_geom_cache(path):
    """geom_cache.pkl 로드. 경로 없거나 파일 없으면 None(→ 호출부가 즉시 파싱으로 폴백)."""
    import pickle
    if path and os.path.exists(path):
        with open(path, "rb") as f:
            return pickle.load(f)
    return None


def cached_geometry(path, cache):
    """
    (prot_ca dict, elems list, coords ndarray(N,3)) 반환.
    - cache 에 path 있으면 파싱 없이 캐시에서 복원(빠름).
    - 없으면 parse_complex + ligand_concat 로 즉시 계산(단독 실행/캐시 실패 폴백).
    반환형은 parse_complex()+ligand_concat() 조합과 동일하게 맞춰 호출부 수정 최소화.
    """
    import numpy as np
    if cache is not None and path in cache:
        e = cache[path]
        prot_ca = {tuple(k): [float(v[0]), float(v[1]), float(v[2])]
                   for k, v in zip(e["pk"], e["px"])}
        return prot_ca, list(e["le"]), np.asarray(e["lx"])
    ca, ligs = parse_complex(path)
    elems, coords = ligand_concat(ligs)
    return ca, elems, (np.asarray(coords) if coords else np.zeros((0, 3)))


def cached_copy_centroids(path, cache):
    """구조의 organic 리간드 copy '각각'의 무게중심 리스트(원본 프레임; 이온 제외).
    - 캐시에 'lc' 있으면 재파싱 없이 사용, 없으면 parse_complex + copy_centroids 폴백.
    - 다중 copy면 리스트 길이 = copy 수(스텝2가 copy별 포켓 후보점으로 사용)."""
    if cache is not None and path in cache and "lc" in cache[path]:
        return [[float(x) for x in c] for c in cache[path]["lc"]]
    _, ligs = parse_complex(path)
    return copy_centroids(ligs)


# 이 파일은 import 전용 공용 라이브러리입니다 (단독 실행 X).
# 스텝들이 `import complex_io` 로 사용 → 스텝들과 같은 폴더에 두어야 함.
