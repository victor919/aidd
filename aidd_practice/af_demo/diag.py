#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 文章关联：AI制药细节篇1_靶点与口袋_公众号文章.html
# 脚本职责：开发诊断：编号、口袋参数、Kabsch 与独立交叉验证
# 使用注解：保留排错过程；不是正文推荐入口，不需要为了阅读文章依次执行。
# 完整对应关系与运行顺序：../README.md
"""诊断：① 用序列比对实测 PDB<->UniProt 编号偏移 ② 扫埋藏度阈值找可用参数"""
import os
import sys
import json

import numpy as np
from scipy.spatial import cKDTree
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.abspath(os.path.join(HERE, "..", ".."))
MVP = os.path.join(WORK, "aidd-mvp")
sys.path.insert(0, os.path.join(MVP, "scripts"))
sys.path.insert(0, os.path.join(MVP, "vendor"))
from _common import raw  # noqa: E402
from e1_af_vs_crystal import parse_pdb  # noqa: E402

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}

c_atoms, c_res = parse_pdb(raw("1M17.pdb"))
a_atoms, a_res = parse_pdb(os.path.join(HERE, "struct", "AF-P00533-F1.pdb"))

print("=== 1M17 晶体结构 ===")
chains = sorted({l[21] for l in open(raw("1M17.pdb")) if l.startswith("ATOM")})
print("  chains:", chains)
nums = sorted(n for n, r in c_res.items() if not r["hetero"])
print("  残基号范围:", nums[0], "-", nums[-1], " 共", len(nums))
resns = sorted({r["resn"] for n, r in c_res.items() if not r["hetero"]})
print("  残基类型数:", len(resns))

print("=== AF-P00533-F1 ===")
anums = sorted(a_res)
print("  残基号范围:", anums[0], "-", anums[-1], " 共", len(anums))

# 扫偏移：offset 使 AF[pdb+offset] 的残基名 == 晶体[pdb] 的残基名
print("=== 偏移扫描（AF残基号 = PDB残基号 + offset） ===")
best = []
for off in range(-40, 41):
    match = tot = 0
    for n in nums:
        m = n + off
        if m in a_res:
            tot += 1
            if a_res[m]["resn"] == c_res[n]["resn"]:
                match += 1
    if tot > 50:
        best.append((match / tot, off, match, tot))
best.sort(reverse=True)
for rate, off, m, t in best[:5]:
    print("  offset=%+d  匹配 %d/%d = %.4f" % (off, m, t, rate))

# 链信息（用于后续只取 chain A）
print("=== 链分布 ===")
for ch in chains:
    sub = [n for n in nums if n in c_res]
    print("  chain", ch, "残留基数(合并前难以区分，需重解析)")

# 埋藏度阈值扫描（1M17，全部蛋白原子，删配体）
def buriedness(coords, spacing=1.0, occ_r=1.8, ray_len=12):
    lo = coords.min(axis=0) - 6.0
    hi = coords.max(axis=0) + 6.0
    nx, ny, nz = [int(x) for x in (np.ceil((hi - lo) / spacing) + 1)]
    gx, gy, gz = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    pts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], 1) * spacing + lo
    tree = cKDTree(coords)
    d, _ = tree.query(pts, k=1)
    occ = (d <= occ_r).reshape(nx, ny, nz)
    AXES = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, -1, 0), (1, 0, 1),
            (1, 0, -1), (0, 1, 1), (0, 1, -1), (1, 1, 1), (1, 1, -1),
            (1, -1, 1), (1, -1, -1)]
    steps = int(ray_len / spacing)
    buried = np.zeros(occ.shape, dtype=np.int16)
    for ax in AXES:
        hp = np.zeros(occ.shape, bool); hn = np.zeros(occ.shape, bool)
        cp = occ.astype(np.float32); cn = occ.astype(np.float32)
        for _ in range(steps):
            cp = ndimage.shift(cp, shift=(-ax[0], -ax[1], -ax[2]), order=0, cval=0.0)
            cn = ndimage.shift(cn, shift=(ax[0], ax[1], ax[2]), order=0, cval=0.0)
            hp |= (cp > 0.5); hn |= (cn > 0.5)
        buried += (hp & hn).astype(np.int16)
    return occ, buried, lo, spacing

heavy = np.array([a["xyz"] for a in c_atoms
                  if not a["atomname"].startswith("H")
                  and a["resn"] not in ("AQ4", "HOH", "WAT", "SO4", "PO4", "GOL", "EDO", "PEG")])
occ, buried, lo, spacing = buriedness(heavy)
free = ~occ
print("=== 埋藏度分布（1M17 空格点） ===")
vals = buried[free]
for t in range(0, 14):
    print("  buriedness>=%2d : %8d 点 (%.2f%%)" % (t, int((vals >= t).sum()), 100 * (vals >= t).mean()))
print("=== 阈值扫描：聚类后的口袋数 ===")
for t in (7, 8, 9, 10, 11, 12, 13):
    mask = free & (buried >= t)
    lab, n = ndimage.label(mask, structure=ndimage.generate_binary_structure(3, 1))
    sizes = ndimage.sum(mask, lab, range(1, n + 1))
    sizes = np.array(sizes)
    big = sizes[sizes >= 30]
    print("  t=%2d 簇数=%4d  >=30埃^3的簇=%3d  最大簇=%7d 埃^3" %
          (t, n, len(big), int(big.max()) if len(big) else 0))
