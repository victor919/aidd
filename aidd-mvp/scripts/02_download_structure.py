"""02 · 下载实验结构并做结构质检

输入：config.structure.pdb_id
输出：data/raw/1M17.pdb + data/processed/02_structure_qc.json

要点：有可靠实验结构时，优先用实验结构，不要一上来就 AlphaFold。
"""
import os
import re
import requests

from _common import cfg, ensure_dirs, raw, processed, log, banner, save_json, exists

ensure_dirs()
banner("02 · Structure：下载实验结构 + 结构质检")

pdb_id = cfg["structure"]["pdb_id"]
lig = cfg["structure"]["ligand_resname"]
url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
dst = raw(f"{pdb_id}.pdb")

if exists(dst):
    log(f"复用已有结构 {os.path.basename(dst)}（{os.path.getsize(dst)/1024:.0f} KB）")
else:
    log(f"GET {url}")
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    with open(dst, "wb") as f:
        f.write(r.content)
    log(f"downloaded {pdb_id} -> data/raw/{pdb_id}.pdb ({len(r.content)/1024:.0f} KB)")

# ---------- 结构质检（纯文本解析，不依赖 Biopython） ----------
txt = open(dst, encoding="utf-8", errors="ignore").read()
lines = txt.splitlines()

res_line = next((l for l in lines if l.startswith("REMARK   2 RESOLUTION")), "")
resolution = None
if res_line:
    m = re.search(r"([0-9]+\.[0-9]+)", res_line[21:])
    resolution = float(m.group(1)) if m else None

atoms = [l for l in lines if l.startswith("ATOM")]
hetatms = [l for l in lines if l.startswith("HETATM")]
chains = sorted({l[21] for l in atoms})
het_res = sorted({l[17:20].strip() for l in hetatms})

# 配体原子数 / 是否有构象 altloc
lig_atoms = [l for l in hetatms if l[17:20].strip() == lig]

# 缺失残基提示
missing = [l[10:].strip() for l in lines if l.startswith("REMARK 465")][:12]

print(f"\n  分辨率          : {resolution} A" if resolution else "\n  分辨率: 未知")
print(f"  蛋白 ATOM 行数  : {len(atoms)}")
print(f"  HETATM 行数     : {len(hetatms)}")
print(f"  蛋白链          : {chains}")
print(f"  非标准残基/配体 : {het_res}")
print(f"  配体 {lig} 原子数  : {len(lig_atoms)}")
if missing:
    print(f"  缺失残基提示    : {len(missing)} 条（示例：{missing[:3]}）")

# ---------- 判定 ----------
qc = {
    "pdb_id": pdb_id,
    "resolution_A": resolution,
    "n_protein_atoms": len(atoms),
    "n_hetatm": len(hetatms),
    "chains": chains,
    "het_residues": het_res,
    "ligand_resname": lig,
    "n_ligand_atoms": len(lig_atoms),
    "has_cognate_ligand": bool(lig_atoms),
}
verdict = []
if resolution is None or resolution > 3.0:
    verdict.append("分辨率偏粗，对接结论要谨慎")
if not lig_atoms:
    verdict.append(f"未找到共晶配体 {lig}，口袋中心需要另找依据")
qc["qc_notes"] = verdict or ["结构与配体齐备，可用于口袋定义"]

save_json(qc, processed("02_structure_qc.json"))
banner(f"02 done · {pdb_id} 分辨率 {resolution} A，含共晶配体 {lig}")
