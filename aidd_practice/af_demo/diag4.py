#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 文章关联：AI制药细节篇1_靶点与口袋_公众号文章.html
# 脚本职责：开发诊断：编号、口袋参数、Kabsch 与独立交叉验证
# 使用注解：保留排错过程；不是正文推荐入口，不需要为了阅读文章依次执行。
# 完整对应关系与运行顺序：../README.md
"""诊断4：用 scipy 独立实现交叉验证 AF vs 晶体的真实偏差"""
import os
import sys

import numpy as np
from scipy.spatial.transform import Rotation

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.abspath(os.path.join(HERE, "..", ".."))
MVP = os.path.join(WORK, "aidd-mvp")
sys.path.insert(0, os.path.join(MVP, "scripts"))
sys.path.insert(0, os.path.join(MVP, "vendor"))
from _common import raw  # noqa: E402
from e1_af_vs_crystal import parse_pdb, fit_transform, OFFSET, CUTOFF, LIGAND_RESN  # noqa: E402


def kabsch_rmsd(P, Q):
    Rm, t = fit_transform(P, Q)
    return float(np.sqrt((((P @ Rm + t) - Q) ** 2).sum(1).mean()))


def scipy_rmsd(P, Q):
    Pc = P - P.mean(0)
    Qc = Q - Q.mean(0)
    rot, rssd = Rotation.align_vectors(Qc, Pc)
    return float(rssd / np.sqrt(len(P)))


c_atoms, c_res = parse_pdb(raw("1M17.pdb"))
a_atoms, a_res = parse_pdb(os.path.join(HERE, "struct", "egfr_af_v6.pdb"))
af_by = {rn + OFFSET: r for rn, r in a_res.items()}
plddt = {a["resnum"]: a["bfac"] for a in a_atoms if a["atomname"] == "CA"}

nums = [n for n in sorted(c_res)
        if n in af_by and c_res[n]["ca"] is not None and not c_res[n]["hetero"]]
P = np.array([af_by[n]["ca"] for n in nums])
Q = np.array([c_res[n]["ca"] for n in nums])
print("共有残基:", len(nums))
print("自有 Kabsch RMSD : %.3f 埃" % kabsch_rmsd(P, Q))
print("scipy align RMSD : %.3f 埃" % scipy_rmsd(P, Q))

# 迭代剔除离群点
idx = np.arange(len(nums))
for it in range(6):
    Rm, t = fit_transform(P[idx], Q[idx])
    dev = np.linalg.norm((P[idx] @ Rm + t) - Q[idx], axis=1)
    r = float(np.sqrt((dev ** 2).mean()))
    keep = max(int(len(idx) * 0.8), 40)
    idx = idx[np.argsort(dev)[:keep]]
    print(f"  第{it+1}轮: 保留{len(idx)}个残基  RMSD={r:.3f} 埃")

# 高置信区
hi = [i for i, n in enumerate(nums) if plddt.get(n - OFFSET, 0) >= 70]
print("pLDDT>=70 残基数:", len(hi))
if len(hi) > 20:
    print("  高置信区 RMSD = %.3f 埃" % kabsch_rmsd(P[hi], Q[hi]))


def dm(X):
    return np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))


D1, D2 = dm(P[hi]), dm(Q[hi])
iu = np.triu_indices(len(hi), 1)
r = np.corrcoef(D1[iu], D2[iu])[0, 1]
rmse = float(np.sqrt(((D1[iu] - D2[iu]) ** 2).mean()))
print(f"  高置信区内部距离 Pearson r = {r:.4f}  RMSE = {rmse:.2f} 埃")

rg_af = float(np.sqrt(((P[hi] - P[hi].mean(0)) ** 2).sum(1).mean()))
rg_cry = float(np.sqrt(((Q[hi] - Q[hi].mean(0)) ** 2).sum(1).mean()))
print(f"  Rg: AF {rg_af:.1f} / 晶体 {rg_cry:.1f} 埃")