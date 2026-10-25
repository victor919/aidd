#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 文章关联：AI制药细节篇1_靶点与口袋_公众号文章.html
# 脚本职责：正文第六节入口：下载、裁剪、找口袋、交互展示
# 使用注解：独立小闭环，输出 struct/、results/e6_end2end.json 和 viz/e6_end2end_view.html；不自动产出全部正文图。
# 完整对应关系与运行顺序：../README.md
"""
E6｜端到端闭环：下载 → 裁剪 → 找口袋 → 出图

单文件可独立运行（只依赖 numpy / scipy / requests / py3Dmol）：

    python e6_end2end.py

四步全部给出真实产物，不是示意：
  step1 下载 AlphaFold 模型（按 UniProt）+ 实验结构（按 PDB ID），落盘到 struct/
  step2 读 pLDDT、按结构域区间裁剪，写出裁剪后的 PDB
  step3 几何口袋检测（LIGSITE 式埋藏度 + 连通聚类），输出 top3 口袋的中心/体积/埋藏度
  step4 用 py3Dmol 把「蛋白 + 口袋球 + 裁剪前后对照」写成单个 HTML

文章第六节引用的就是本脚本的四段核心代码，变量全部在本文件内定义。
"""
import os
import sys
import json
import argparse
import urllib.request

import numpy as np
from scipy.spatial import cKDTree
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(WORK, "aidd-mvp", "scripts"))
sys.path.insert(0, os.path.join(WORK, "aidd-mvp", "vendor"))

STRUCT = os.path.join(HERE, "struct")
RESULTS = os.path.join(HERE, "results")
FIG = os.path.join(WORK, "AI制药细节篇1_配图")
os.makedirs(STRUCT, exist_ok=True)
os.makedirs(RESULTS, exist_ok=True)

UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/120.0 Safari/537.36")

# ---------------- 通用：PDB 解析 ----------------

def parse_pdb(path):
    """返回 (atoms, residues)。atoms: 每行一个原子；residues: 残基号 -> {resn, atoms, ca}"""
    residues, atoms = {}, []
    for line in open(path):
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        if line[16] not in (" ", "A"):
            continue
        resname = line[17:20].strip()
        try:
            resnum = int(line[22:26])
        except ValueError:
            continue
        atomname = line[12:16].strip()
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        bfac = float(line[60:66]) if line[60:66].strip() else 0.0
        atoms.append(dict(resnum=resnum, resn=resname, atomname=atomname,
                          xyz=xyz, bfac=bfac, line=line))
        r = residues.setdefault(resnum, dict(resn=resname, atoms={}, hetero=False))
        r["atoms"].setdefault(atomname, xyz)
        if line.startswith("HETATM"):
            r["hetero"] = True
    for r in residues.values():
        r["ca"] = r["atoms"].get("CA")
    return atoms, residues


def band(plddt):
    """pLDDT 分档（AlphaFold DB 官方四档阈值）"""
    if plddt >= 90:
        return "very_high"
    if plddt >= 70:
        return "confident"
    if plddt >= 50:
        return "low"
    return "very_low"


# ---------------- step 1：下载 ----------------

def fetch_alphafold(uniprot, out_path):
    if os.path.exists(out_path) and os.path.getsize(out_path) > 10000:
        print("  [缓存] %s (%.0f KB)" % (os.path.basename(out_path),
                                         os.path.getsize(out_path) / 1024))
        return out_path
    api = "https://alphafold.ebi.ac.uk/api/prediction/%s" % uniprot
    req = urllib.request.Request(api, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=60) as r:
        data = json.loads(r.read().decode())
    entry = data[0] if isinstance(data, list) else data
    url = entry["pdbUrl"]
    req = urllib.request.Request(url, headers={"User-Agent": UA,
                                               "Referer": "https://alphafold.ebi.ac.uk/"})
    with urllib.request.urlopen(req, timeout=120) as r:
        txt = r.read().decode()
    open(out_path, "w").write(txt)
    print("  [下载] AlphaFold %s → %s (%.0f KB, %s)"
          % (uniprot, os.path.basename(out_path), len(txt) / 1024,
             entry.get("latestVersion", "?")))
    return out_path


def fetch_rcsb(pdb_id, out_path):
    if os.path.exists(out_path) and os.path.getsize(out_path) > 10000:
        print("  [缓存] %s (%.0f KB)" % (os.path.basename(out_path),
                                         os.path.getsize(out_path) / 1024))
        return out_path
    url = "https://files.rcsb.org/download/%s.pdb" % pdb_id
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=60) as r:
        txt = r.read().decode()
    open(out_path, "w").write(txt)
    print("  [下载] RCSB %s → %s (%.0f KB)" % (pdb_id, os.path.basename(out_path),
                                               len(txt) / 1024))
    return out_path


# ---------------- step 2：读 pLDDT + 裁剪 ----------------

def read_plddt_and_crop(af_path, uniprot_lo, uniprot_hi, out_pdb):
    atoms, residues = parse_pdb(af_path)
    plddt = {a["resnum"]: a["bfac"] for a in atoms if a["atomname"] == "CA"}

    counts = {}
    for v in plddt.values():
        counts[band(v)] = counts.get(band(v), 0) + 1
    total = len(plddt)
    print("  pLDDT 分档（共 %d 个残基）:" % total)
    for b in ("very_high", "confident", "low", "very_low"):
        n = counts.get(b, 0)
        print("    %-11s %4d 个  %5.1f%%" % (b, n, 100.0 * n / total))

    keep = {n for n in residues if uniprot_lo <= n <= uniprot_hi}
    wrote = 0
    with open(out_pdb, "w") as f:
        f.write("REMARK   cropped from %s, UniProt %d-%d\n"
                % (os.path.basename(af_path), uniprot_lo, uniprot_hi))
        for a in atoms:
            if a["resnum"] in keep:
                f.write(a["line"] if a["line"].endswith("\n") else a["line"] + "\n")
                wrote += 1
        f.write("END\n")
    kept_plddt = [plddt[n] for n in keep if n in plddt]
    print("  裁剪 UniProt %d-%d：保留 %d 个残基 / %d 个原子"
          % (uniprot_lo, uniprot_hi, len(keep), wrote))
    if kept_plddt:
        print("  裁剪区 pLDDT：均值 %.2f，中位数 %.2f，最低 %.2f"
              % (np.mean(kept_plddt), np.median(kept_plddt), np.min(kept_plddt)))
    return dict(band_counts=counts, n_total=total, n_kept_residues=len(keep),
                n_kept_atoms=wrote,
                kept_plddt_mean=round(float(np.mean(kept_plddt)), 2),
                kept_plddt_median=round(float(np.median(kept_plddt)), 2),
                kept_plddt_min=round(float(np.min(kept_plddt)), 2))


# ---------------- step 3：几何口袋检测 ----------------

SPACING, OCC_R, RAY_LEN, BURY_MIN, MIN_CLUSTER = 1.0, 1.8, 12, 9, 50
SOLVENTS = {"HOH", "WAT", "SO4", "PO4", "GOL", "EDO", "PEG", "CL", "NA", "MG", "ZN"}
AXES = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, -1, 0), (1, 0, 1),
        (1, 0, -1), (0, 1, 1), (0, 1, -1), (1, 1, 1), (1, 1, -1),
        (1, -1, 1), (1, -1, -1)]


def detect_pockets(pdb_path, exclude_resn=()):
    atoms, _ = parse_pdb(pdb_path)
    heavy = np.array([a["xyz"] for a in atoms
                      if not a["atomname"].startswith("H")
                      and a["resn"] not in SOLVENTS
                      and a["resn"] not in exclude_resn])

    lo = heavy.min(axis=0) - 6.0
    hi = heavy.max(axis=0) + 6.0
    nx, ny, nz = [int(x) for x in (np.ceil((hi - lo) / SPACING) + 1)]
    gx, gy, gz = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    pts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], 1) * SPACING + lo
    d, _ = cKDTree(heavy).query(pts, k=1)
    occ = (d <= OCC_R).reshape(nx, ny, nz)

    buried = np.zeros(occ.shape, dtype=np.int16)
    for ax in AXES:
        hp = np.zeros(occ.shape, bool)
        hn = np.zeros(occ.shape, bool)
        cp = occ.astype(np.float32)
        cn = occ.astype(np.float32)
        for _ in range(int(RAY_LEN / SPACING)):
            cp = ndimage.shift(cp, shift=(-ax[0], -ax[1], -ax[2]), order=0, cval=0.0)
            cn = ndimage.shift(cn, shift=(ax[0], ax[1], ax[2]), order=0, cval=0.0)
            hp |= (cp > 0.5)
            hn |= (cn > 0.5)
        buried += (hp & hn).astype(np.int16)

    mask = (~occ) & (buried >= BURY_MIN)
    lab, nlab = ndimage.label(mask, structure=ndimage.generate_binary_structure(3, 1))
    pockets = []
    for i in range(1, nlab + 1):
        sel = (lab == i)
        vol = int(sel.sum())
        if vol < MIN_CLUSTER:
            continue
        idx = np.argwhere(sel)
        pockets.append(dict(volume_A3=vol,
                            mean_buriedness=round(float(buried[sel].mean()), 2),
                            center=[round(float(x), 2) for x in idx.mean(axis=0) * SPACING + lo]))
    pockets.sort(key=lambda p: -p["volume_A3"] * p["mean_buriedness"])
    return pockets, heavy


# ---------------- step 4：出图 ----------------

def write_viewer(af_full, af_crop, pockets, out_html):
    import py3Dmol
    v = py3Dmol.view(width=880, height=560)
    v.addModel(open(af_full).read(), "pdb")
    v.setStyle({}, {"cartoon": {"color": "#d8dee9"}})
    v.addModel(open(af_crop).read(), "pdb")
    v.setStyle({"model": -1}, {"cartoon": {"color": "#0FA87E"}})
    for i, p in enumerate(pockets[:3]):
        c = p["center"]
        v.addSphere({"center": {"x": c[0], "y": c[1], "z": c[2]},
                     "radius": 4.0,
                     "color": ["#E53E3E", "#F5A623", "#4A7CFE"][i],
                     "alpha": 0.45})
    v.zoomTo()
    html = v.write_html()
    # 让文件可离线转发：把 CDN 换成本地相对路径不可行，这里保留 CDN 但补一句提示
    open(out_html, "w").write(html)
    return out_html


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--uniprot", default="P00533")
    ap.add_argument("--pdb", default="1M17")
    ap.add_argument("--ligand", default="AQ4")
    ap.add_argument("--lo", type=int, default=696, help="目标结构域 UniProt 起点")
    ap.add_argument("--hi", type=int, default=1019, help="目标结构域 UniProt 终点")
    args = ap.parse_args()

    print("=" * 60)
    print("step 1  下载结构")
    print("=" * 60)
    af_full = fetch_alphafold(args.uniprot, os.path.join(STRUCT, "egfr_af_v6.pdb"))
    crystal = fetch_rcsb(args.pdb, os.path.join(STRUCT, "%s.pdb" % args.pdb))

    print()
    print("=" * 60)
    print("step 2  读 pLDDT 并按结构域裁剪")
    print("=" * 60)
    af_crop = os.path.join(STRUCT, "e6_af_kinase_domain.pdb")
    s2 = read_plddt_and_crop(af_full, args.lo, args.hi, af_crop)

    print()
    print("=" * 60)
    print("step 3  在裁剪后的结构上找口袋")
    print("=" * 60)
    pockets, heavy = detect_pockets(af_crop)
    print("  蛋白重原子 %d 个，检出 %d 个候选口袋（体积 ≥ %d 立方埃）"
          % (len(heavy), len(pockets), MIN_CLUSTER))
    for i, p in enumerate(pockets[:3], 1):
        print("    #%d  体积 %6d 立方埃   平均埋藏度 %.2f   中心 [%.2f, %.2f, %.2f]"
              % (i, p["volume_A3"], p["mean_buriedness"], *p["center"]))

    print()
    print("=" * 60)
    print("step 4  出图")
    print("=" * 60)
    out_html = os.path.join(HERE, "viz", "e6_end2end_view.html")
    os.makedirs(os.path.dirname(out_html), exist_ok=True)
    write_viewer(af_full, af_crop, pockets, out_html)
    print("  交互视图 → %s (%.0f KB)" % (out_html, os.path.getsize(out_html) / 1024))

    out = dict(step1=dict(af=os.path.basename(af_full), crystal=os.path.basename(crystal)),
               step2=s2,
               step3=dict(n_protein_atoms=int(len(heavy)),
                          n_pockets=len(pockets), top3=pockets[:3]),
               step4=dict(html=os.path.basename(out_html)))
    with open(os.path.join(RESULTS, "e6_end2end.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print()
    print("写入 results/e6_end2end.json")


if __name__ == "__main__":
    main()
