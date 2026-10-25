"""03 · 提取口袋残基 + 残基编号偏移校验

输入：data/raw/1M17.pdb、UniProt P00533 序列（联网获取并缓存）
输出：data/processed/pocket.json、pocket_residues.csv、03_pocket.json

要点（本项目最大的坑）：
  1M17 的残基编号 = UniProt 编号 - 24。结构里的 C773 就是临床说的 C797。
  不做这一步校验，你后面所有"关键残基"结论都会错。
"""
import os
import re
import json
import requests
from collections import Counter

from _common import cfg, ensure_dirs, raw, processed, log, banner, save_json, exists

ensure_dirs()
banner("03 · Pocket：配体周围 5 A 残基 + 编号偏移校验")

from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa

pdb_id = cfg["structure"]["pdb_id"]
lig = cfg["structure"]["ligand_resname"]
radius = float(cfg["structure"]["pocket_radius"])
offset_cfg = int(cfg["structure"]["numbering_offset_to_uniprot"])
uniprot = cfg["target"]["uniprot"]

pdb_path = raw(f"{pdb_id}.pdb")
if not exists(pdb_path):
    raise SystemExit("请先运行 02_download_structure.py")

# ---------- 1. 邻居搜索：配体 radius 内的蛋白残基 ----------
parser = PDBParser(QUIET=True)
st = parser.get_structure("T", pdb_path)[0]      # 只取第一个 model

lig_atoms = []
for res in st.get_residues():
    if res.get_resname().strip() == lig:
        lig_atoms.extend(list(res.get_atoms()))
log(f"配体 {lig} 原子 {len(lig_atoms)} 个")

prot_atoms = [a for a in st.get_atoms()
              if a.get_parent().id[0] == " " and is_aa(a.get_parent().get_resname())]

ns = NeighborSearch(prot_atoms)
contacts = Counter()
for a in lig_atoms:
    for r in ns.search(a.coord, radius, level="R"):
        contacts[(r.get_parent().id, r.id[1], r.get_resname())] += 1

res_list = sorted(contacts.items(), key=lambda kv: -kv[1])
log(f"{radius} A 内残基 {len(res_list)} 个")

print(f"\n{'链':<3}{'编号(PDB)':<10}{'残基':<6}{'接触原子对':<10}{'UniProt 编号'}")
print("-" * 48)
for (ch, num, rn), c in res_list:
    print(f"{ch:<3}{num:<10}{rn:<6}{c:<10}{num + offset_cfg}")

# ---------- 2. 编号偏移校验：拿结构序列去比对 UniProt ----------
fasta_path = raw(f"{uniprot}.fasta")
if not exists(fasta_path):
    log(f"GET UniProt {uniprot} 序列")
    r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot}.fasta", timeout=60)
    r.raise_for_status()
    with open(fasta_path, "w", encoding="utf-8") as f:
        f.write(r.text)
seq = "".join(l.strip() for l in open(fasta_path) if not l.startswith(">"))
log(f"UniProt {uniprot} 序列长度 {len(seq)}")

# 结构里每个链的观测残基
three2one = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}
best = None
for ch in st:
    obs = {}
    for r in ch:
        if r.id[0] == " " and r.get_resname() in three2one:
            obs[r.id[1]] = three2one[r.get_resname()]
    if len(obs) < 50:
        continue
    # 枚举偏移量，找匹配率最高的
    for off in range(-40, 41):
        hit = sum(1 for n, aa in obs.items() if 0 <= n + off - 1 < len(seq)
                  and seq[n + off - 1] == aa)
        rate = hit / len(obs)
        if best is None or rate > best[2]:
            best = (ch.id, off, rate, hit, len(obs))

ch_id, off_detected, rate, hit, n_obs = best
print(f"\n编号偏移校验（链 {ch_id}）：")
print(f"  观测残基 {n_obs} 个，最佳偏移 = {off_detected}，序列匹配 {hit}/{n_obs} = {rate:.3%}")
print(f"  config 里写的偏移 = {offset_cfg}  ->  {'一致' if off_detected == offset_cfg else '不一致，请改 config'}")

# ---------- 3. 用 UniProt 编号复核经典残基 ----------
canonical = {718: "L718", 719: "G719", 745: "K745 催化赖氨酸", 790: "T790 gatekeeper",
             793: "M793 hinge", 797: "C797 共价残基", 844: "L844",
             854: "T854", 855: "D855 DFG-Asp"}
pdb_nums = {n for (_, n, _) in contacts.keys()}
up_nums = {n + off_detected for n in pdb_nums}
print(f"\n经典 EGFR 激酶位点是否在 {radius} A 口袋内（按 UniProt 编号）：")
hits = {}
for num, label in canonical.items():
    ok = num in up_nums
    hits[label] = ok
    print(f"  {label:<20s} {'YES' if ok else 'no'}")

save_json({
    "pdb_id": pdb_id, "ligand": lig, "radius": radius,
    "pocket_residues": [
        {"chain": c, "pdb_num": int(n), "resname": rn,
         "uniProt_num": int(n + off_detected), "contacts": int(k)}
        for ((c, n, rn), k) in res_list],
    "numbering": {"chain": ch_id, "detected_offset": int(off_detected),
                  "config_offset": offset_cfg,
                  "matched": f"{hit}/{n_obs}", "match_rate": round(rate, 4)},
    "canonical_sites_in_pocket": hits,
}, processed("pocket.json"))

import pandas as pd
pd.DataFrame([{"chain": c, "pdb_num": n, "resname": rn,
               "uniprot_num": n + off_detected, "contacts": k}
              for ((c, n, rn), k) in res_list]).to_csv(
    processed("pocket_residues.csv"), index=False)

banner(f"03 done · 口袋 {len(res_list)} 残基，编号偏移确认 = {off_detected}")
