#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 文章关联：AI制药细节篇1_靶点与口袋_公众号文章.html
# 脚本职责：开发诊断：编号、口袋参数、Kabsch 与独立交叉验证
# 使用注解：保留排错过程；不是正文推荐入口，不需要为了阅读文章依次执行。
# 完整对应关系与运行顺序：../README.md
"""诊断2：① E1 的 common 残基为何为空 ② 扫口袋检测参数，找能区分的参数组合"""
import os
import sys

import numpy as np
from scipy.spatial import cKDTree
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.abspath(os.path.join(HERE, "..", ".."))
MVP = os.path.join(WORK, "aidd-mvp")
sys.path.insert(0, os.path.join(MVP, "scripts"))
sys.path.insert(0, os.path.join(MVP, "vendor"))
from _common import raw  # noqa: E402
from e1_af_vs_crystal import parse_pdb, OFFSET, CUTOFF, LIGAND_RESN  # noqa: E402

print("OFFSET =", OFFSET)
c_atoms, c_res = parse_pdb(raw("1M17.pdb"))
a_atoms, a_res = parse_pdb(os.path.join(HERE, "struct", "AF-P00533-F1.pdb"))
print("c_res 个数:", len(c_res), " a_res 个数:", len(a_res))
print("c_res 有 CA 的:", sum(1 for r in c_res.values() if r["ca"] is not None))
print("a_res 有 CA 的:", sum(1 for r in a_res.values() if r["ca"] is not None))
print("c_res 非 hetero:", sum(1 for r in c_res.values() if not r["hetero"]))
print("a_res 非 hetero:", sum(1 for r in a_res.values() if not r["hetero"]))
cn = sorted(c_res)[:5], sorted(c_res)[-5:]
an = sorted(a_res)[:5], sorted(a_res)[-5:]
print("c_res 号样本:", cn)
print("a_res 号样本:", an)
af_keys = {rn + OFFSET for rn in a_res}
common = [n for n in sorted(c_res) if n in af_keys and not c_res[n]["hetero"]]
print("common 个数:", len(common), "样本:", common[:5], common[-5:])
print("样例比对: PDB", common[:3], "-> AF", [n - OFFSET for n in common[:3]])
for n in common[:3]:
    print("   PDB %d %s   AF %d %s" % (n, c_res[n]["resn"], n - OFFSET,
                                       a_res[n - OFFSET]["resn"]))

# ---- 参数扫描
lig = [a for a in c_atoms if a["resn"] == LIGAND_RESN]
lig_xyz = np.array([a["xyz"] for a in lig])
lig_center = lig_xyz.mean(axis=0)
true_resnums = set()
for rnum, r in c_res.items():
    if r["hetero"] or r["ca"] is None:
        continue
    if any(np.linalg.norm(v - lig_center) <= CUTOFF for v in r["atoms"].values()):
        true_resnums.add(rnum)
xyz = [c_res[n]["ca"] for n in true_resnums]
true_center = np.array(xyz).mean(axis=0)
print("真位点残基数:", len(true_resnums), "质心:", np.round(true_center, 2))

heavy = np.array([a["xyz"] for a in c_atoms
                  if not a["atomname"].startswith("H")
                  and a["resn"] not in (LIGAND_RESN, "HOH", "WAT", "SO4", "PO4", "GOL", "EDO")])
AXES = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, -1, 0), (1, 0, 1),
        (1, 0, -1), (0, 1, 1), (0, 1, -1), (1, 1, 1), (1, 1, -1),
        (1, -1, 1), (1, -1, -1)]


def sweep(occ_r, ray_len, bury_min, spacing=1.0, min_cluster=50):
    lo = heavy.min(axis=0) - 6.0
    hi = heavy.max(axis=0) + 6.0
    nx, ny, nz = [int(x) for x in (np.ceil((hi - lo) / spacing) + 1)]
    gx, gy, gz = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    pts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], 1) * spacing + lo
    d, _ = cKDTree(heavy).query(pts, k=1)
    occ = (d <= occ_r).reshape(nx, ny, nz)
    steps = int(ray_len / spacing)
    buried = np.zeros(occ.shape, dtype=np.int16)
    for ax in AXES:
        hp = np.zeros(occ.shape, bool); hn = np.zeros(occ.shape, bool)
        cp = occ.astype(np.float32); cn_ = occ.astype(np.float32)
        for _ in range(steps):
            cp = ndimage.shift(cp, shift=(-ax[0], -ax[1], -ax[2]), order=0, cval=0.0)
            cn_ = ndimage.shift(cn_, shift=(ax[0], ax[1], ax[2]), order=0, cval=0.0)
            hp |= (cp > 0.5); hn |= (cn_ > 0.5)
        buried += (hp & hn).astype(np.int16)
    mask = (~occ) & (buried >= bury_min)
    lab, nlab = ndimage.label(mask, structure=ndimage.generate_binary_structure(3, 1))
    pk = []
    for i in range(1, nlab + 1):
        sel = (lab == i)
        vol = int(sel.sum())
        if vol < min_cluster:
            continue
        cen = np.argwhere(sel).mean(axis=0) * spacing + lo
        mb = float(buried[sel].mean())
        pk.append((vol * mb, vol, mb, float(np.linalg.norm(cen - true_center)), cen))
    pk.sort(key=lambda x: -x[0])
    if not pk:
        return 0, None, None, None
    dccs = [p[3] for p in pk]
    best_i = int(np.argmin(dccs))
    return len(pk), pk[best_i][3], best_i + 1, pk[0][1]


print("\n=== 参数扫描（occ_r / ray_len / bury_min） ===")
print("%-6s %-8s %-9s | %-5s %-9s %-6s %s" % ("occ_r", "ray_len", "bury_min", "簇数", "最优DCC", "排名", "top1体积"))
for occ_r in (1.8, 2.4, 3.0):
    for ray_len in (8, 10, 12):
        for bury_min in (10, 12, 13):
            n, dcc, rank, v1 = sweep(occ_r, ray_len, bury_min)
            print("%-6.1f %-8d %-9d | %-5s %-9s %-6s %s" % (
                occ_r, ray_len, bury_min, n,
                "%.2f" % dcc if dcc is not None else "-", rank if rank else "-", v1))
