#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 文章关联：AI制药细节篇1_靶点与口袋_公众号文章.html
# 脚本职责：开发诊断：编号、口袋参数、Kabsch 与独立交叉验证
# 使用注解：保留排错过程；不是正文推荐入口，不需要为了阅读文章依次执行。
# 完整对应关系与运行顺序：../README.md
"""诊断3：① 用合成数据验证 Kabsch ② 查真实 RMSD 到底是多少"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.abspath(os.path.join(HERE, "..", ".."))
MVP = os.path.join(WORK, "aidd-mvp")
sys.path.insert(0, os.path.join(MVP, "scripts"))
sys.path.insert(0, os.path.join(MVP, "vendor"))
from _common import raw  # noqa: E402
from e1_af_vs_crystal import parse_pdb, fit_transform, OFFSET, CUTOFF, LIGAND_RESN  # noqa: E402

rng = np.random.default_rng(0)
# 合成测试：Q 是随机点云，P 是 Q 经过已知旋转+平移+噪声
Q = rng.normal(size=(100, 3)) * 10
ang = 0.7
Rz = np.array([[np.cos(ang), -np.sin(ang), 0],
               [np.sin(ang), np.cos(ang), 0], [0, 0, 1]])
P = Q @ Rz.T + np.array([5.0, -3.0, 2.0]) + rng.normal(size=(100, 3)) * 0.2
R, t = fit_transform(P, Q)
aligned = P @ R + t
print("[合成] RMSD after fit = %.4f 埃 (期望 ~0.2)" %
      np.sqrt(((aligned - Q) ** 2).sum(1).mean()))
print("[合成] R 是否接近 Rz^T?", np.allclose(R, Rz.T, atol=1e-2))

# 真实数据
c_atoms, c_res = parse_pdb(raw("1M17.pdb"))
a_atoms, a_res = parse_pdb(os.path.join(HERE, "struct", "egfr_af_v6.pdb"))
af_by_uniprot = {rn + OFFSET: r for rn, r in a_res.items()}
common = [n for n in sorted(c_res)
          if n in af_by_uniprot and c_res[n]["ca"] is not None
          and af_by_uniprot[n]["ca"] is not None and not c_res[n]["hetero"]]
P = np.array([af_by_uniprot[n]["ca"] for n in common])
Q = np.array([c_res[n]["ca"] for n in common])
print("[真实] common =", len(common))
print("[真实] 未对齐直接算 RMSD = %.2f 埃" % np.sqrt(((P - Q) ** 2).sum(1).mean()))
R, t = fit_transform(P, Q)
print("[真实] 对齐后 RMSD = %.2f 埃" % np.sqrt((((P @ R + t) - Q) ** 2).sum(1).mean()))

# 内部距离矩阵对比（与叠合无关，最能说明"折叠对不对"）
def dm(X):
    return np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))
D1, D2 = dm(P), dm(Q)
iu = np.triu_indices(len(common), 1)
print("[真实] 内部距离矩阵 RMSE(Cα) = %.3f 埃" %
      np.sqrt(((D1[iu] - D2[iu]) ** 2).mean()))
print("[真实] 内部距离 Pearson r = %.4f" % np.corrcoef(D1[iu], D2[iu])[0, 1])
rg_af = np.sqrt(((P - P.mean(0)) ** 2).sum(1).mean())
rg_cry = np.sqrt(((Q - Q.mean(0)) ** 2).sum(1).mean())
print("[真实] 半径回转 Rg: AF %.1f 埃 / 晶体 %.1f 埃" % (rg_af, rg_cry))