# 文章关联：AI制药实践篇B_实操篇_公众号文章.html
# 脚本职责：早期实验 1：靶点排序、1M17 与口袋提取
# 使用注解：写入 data/egfr_mvp/；属于正式 aidd-mvp 工程之前的试跑口径。
# 完整对应关系与运行顺序：README.md
"""EGFR MVP Step1: target ranking + 1M17 structure + pocket extraction
所有输出均为本机实跑结果。
"""
import os, sys, json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "vendor"))

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(BASE, "data", "egfr_mvp")
os.makedirs(DATA, exist_ok=True)

print("=" * 70)
print("[1] Target ranking (multi-evidence scoring)")
print("=" * 70)

rows = [
    # target, disease, genetics, pathway, druggability, safety(risk)
    ("EGFR",      0.95, 0.82, 0.90, 0.95, 0.30),
    ("Target_B",  0.70, 0.91, 0.65, 0.50, 0.45),
    ("Target_C",  0.80, 0.55, 0.82, 0.35, 0.25),
]
df = pd.DataFrame(rows, columns=["target", "disease", "genetics", "pathway", "druggability", "safety"])
df["score"] = (0.30 * df["disease"] + 0.25 * df["genetics"] + 0.20 * df["pathway"]
               + 0.25 * df["druggability"] - 0.20 * df["safety"])
df = df.sort_values("score", ascending=False).reset_index(drop=True)
print(df[["target", "score"]].to_string(index=False))
df.to_csv(os.path.join(DATA, "target_rank.csv"), index=False)

# 敏感性检查：换权重，排名是否稳定
print("\n-- weight sensitivity --")
W = {
    "baseline":      (0.30, 0.25, 0.20, 0.25, 0.20),
    "safety-first":  (0.25, 0.20, 0.15, 0.20, 0.40),
    "genetics-first":(0.20, 0.40, 0.15, 0.25, 0.20),
}
for name, (a, b, c, d, e) in W.items():
    s = a*df["disease"] + b*df["genetics"] + c*df["pathway"] + d*df["druggability"] - e*df["safety"]
    order = df.assign(s=s).sort_values("s", ascending=False)["target"].tolist()
    print(f"{name:15s} -> {' > '.join(order)}   scores={[round(x,3) for x in df.assign(s=s)['s'].tolist()]}")

print("\n" + "=" * 70)
print("[2] Download 1M17 and inspect")
print("=" * 70)
import requests
pdb_id = "1M17"
pdb_path = os.path.join(DATA, f"{pdb_id}.pdb")
if not os.path.exists(pdb_path):
    r = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=60)
    r.raise_for_status()
    open(pdb_path, "wb").write(r.content)
print("file:", pdb_path, os.path.getsize(pdb_path), "bytes")

for line in open(pdb_path):
    if line.startswith("REMARK   2 RESOLUTION"):
        print(" ", line.strip())

print("\n" + "=" * 70)
print("[3] Pocket extraction with Biopython NeighborSearch")
print("=" * 70)
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa

LIGAND = "AQ4"
RADIUS = 5.0
parser = PDBParser(QUIET=True)
structure = parser.get_structure("EGFR", pdb_path)
model = structure[0]

ligand_atoms, protein_atoms = [], []
for chain in model:
    for residue in chain:
        if residue.get_resname().strip() == LIGAND:
            ligand_atoms.extend(list(residue.get_atoms()))
        elif is_aa(residue, standard=True):
            protein_atoms.extend(list(residue.get_atoms()))

print(f"ligand {LIGAND}: {len(ligand_atoms)} atoms")
print(f"standard protein residues atoms: {len(protein_atoms)}")

ns = NeighborSearch(protein_atoms)
pocket = {}
for atom in ligand_atoms:
    for res in ns.search(atom.coord, RADIUS, level="R"):
        key = (res.parent.id, res.id[1], res.resname)
        pocket[key] = pocket.get(key, 0) + 1

pocket_sorted = sorted(pocket.items(), key=lambda kv: (kv[0][0], kv[0][1]))
print(f"\npocket residues within {RADIUS} A of {LIGAND}: {len(pocket_sorted)}")
print("chain resnum resname  contacts")
for (ch, num, name), cnt in pocket_sorted:
    print(f"  {ch}   {num:<5d} {name:<4s}  {cnt}")

# 更严格的 4.0 A
ns4 = NeighborSearch(protein_atoms)
p4 = set()
for atom in ligand_atoms:
    for res in ns4.search(atom.coord, 4.0, level="R"):
        p4.add((res.parent.id, res.id[1], res.resname))
print(f"\n(4.0 A cutoff -> {len(p4)} residues)")

# 关键残基核对：注意 1M17 编号 = UniProt P00533 编号 - 24（见 egfr_mvp_numbering.py 100% 校验）
OFFSET = 24
canonical = {"L718": 718, "G719": 719, "K745": 745, "T790": 790, "M793": 793,
             "C797": 797, "L844": 844, "T854": 854, "D855": 855}
have_pdb = {key[1] for key, _ in pocket_sorted}
have_up = {key[1] + OFFSET for key, _ in pocket_sorted}
print("\ncanonical EGFR kinase-site residues (UniProt numbering) found in 5A pocket:")
for label, num in canonical.items():
    print(f"  {label}: {'YES' if num in have_up else 'no'}   (1M17 res {num - OFFSET})")

# pocket 几何中心（docking box 中心）
lig_coords = np.array([a.coord for a in ligand_atoms], dtype=float)
print(f"\nligand centroid (docking box center): {np.round(lig_coords.mean(axis=0), 2).tolist()}")
print(f"ligand bbox size (A): {np.round(lig_coords.max(axis=0)-lig_coords.min(axis=0), 2).tolist()}")

json.dump({
    "pdb": pdb_id, "ligand": LIGAND, "radius": RADIUS,
    "n_pocket_residues_5A": len(pocket_sorted),
    "n_pocket_residues_4A": len(p4),
    "residues": [[ch, int(n), rn] for (ch, n, rn), _ in pocket_sorted],
    "box_center": [round(float(x), 2) for x in lig_coords.mean(axis=0)],
}, open(os.path.join(DATA, "pocket.json"), "w"), indent=2)
print("\nsaved ->", os.path.join(DATA, "pocket.json"))
