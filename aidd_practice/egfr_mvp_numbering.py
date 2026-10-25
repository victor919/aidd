# 文章关联：AI制药实践篇B_实操篇_公众号文章.html
# 脚本职责：事实核查：PDB 与 UniProt 残基编号偏移
# 使用注解：依赖 data/egfr_mvp/1M17.pdb；为口袋编号说明提供校验依据。
# 完整对应关系与运行顺序：README.md
"""核实 1M17 结构编号 vs UniProt P00533 编号的偏移关系（硬事实校验）"""
import os, sys, re, json
import requests

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "vendor"))
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

THREE2ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(BASE, "data", "egfr_mvp")

# ---------- 1. UniProt P00533 序列 ----------
fasta_url = "https://rest.uniprot.org/uniprotkb/P00533.fasta"
r = requests.get(fasta_url, timeout=60)
r.raise_for_status()
lines = r.text.strip().split("\n")
header = lines[0]
seq = "".join(lines[1:])
print("UniProt header:", header[:90])
print("P00533 length :", len(seq))

def aa(pos):  # 1-based
    return seq[pos - 1]

print("\n-- UniProt 编号 (临床常用) 位点核对 --")
for p, expect in [(719, "G"), (745, "K"), (766, None), (790, "T"), (773, None),
                  (797, "C"), (820, None), (844, "L"), (854, "T"), (855, "D"),
                  (858, "L"), (861, "L"), (793, "M"), (796, "G"), (762, "E"), (718, "L")]:
    print(f"  {p}: {aa(p)}   (expect {expect})")

# ---------- 2. 1M17 观测残基 ----------
parser = PDBParser(QUIET=True)
st = parser.get_structure("E", os.path.join(DATA, "1M17.pdb"))
obs = {}
for ch in st[0]:
    for res in ch:
        if is_aa(res, standard=True):
            obs[res.id[1]] = THREE2ONE[res.get_resname().strip().upper()]
print(f"\n1M17 chain A observed residues: {len(obs)}, range {min(obs)}-{max(obs)}")

# ---------- 3. 求偏移 ----------
print("\n-- 偏移检验：obs[k] == seq[k + off]? --")
best = None
for off in range(-40, 41):
    match = sum(1 for k, v in obs.items() if 1 <= k + off <= len(seq) and seq[k + off - 1] == v)
    frac = match / len(obs)
    if best is None or frac > best[1]:
        best = (off, frac)
print(f"  best offset = {best[0]}  match = {best[1]*100:.2f}%  ({int(best[1]*len(obs))}/{len(obs)})")

off = best[0]
print(f"\n-- 映射表 (1M17 -> UniProt, offset {off}) --")
print("  1M17  residue  ->  UniProt  residue   (临床常用名)")
pairs = [(694, "L718"), (695, "G719"), (702, "V726"), (721, "K745"), (766, "T790"),
         (769, "M793"), (772, "G796"), (773, "C797"), (820, "L844"), (830, "T854"),
         (831, "D855"), (719, "?"), (738, "E762"), (742, "M766")]
for pdbnum, clinical in pairs:
    up = pdbnum + off
    ok = (1 <= up <= len(seq))
    print(f"  {pdbnum:<5d} {obs.get(pdbnum,'-'):<4s}   ->   {up:<5d} {aa(up) if ok else '-':<4s}    {clinical}"
          f"   {'MATCH' if ok and clinical[1:-len(clinical[1:])+1:] else ''}")

print("\n-- 临床位点反查 --")
for name, uppos in [("G719", 719), ("K745", 745), ("T790", 790), ("C797", 797),
                    ("L844", 844), ("T854", 854), ("D855", 855), ("L858", 858)]:
    pdbnum = uppos - off
    key = f"{aa(uppos)}{uppos}"
    print(f"  {name} (UniProt {uppos} {aa(uppos)})  ->  1M17 residue {pdbnum} = {obs.get(pdbnum, 'not observed')}"
          f"   {'OK' if obs.get(pdbnum) == aa(uppos) else 'MISMATCH/NA'}")

# ---------- 4. 配体身份 ----------
print("\n-- 1M17 配体 --")
for line in open(os.path.join(DATA, "1M17.pdb")):
    if line.startswith("HET   ") or (line.startswith("HETNAM") and "AQ4" in line):
        print("  " + line.strip())
for line in open(os.path.join(DATA, "1M17.pdb")):
    if line.startswith("TITLE"):
        print("  " + line.strip())
