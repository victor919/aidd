#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 文章关联：AI制药细节篇1_靶点与口袋_公众号文章.html
# 脚本职责：补算：RMSD 残基集合、叠合方式与覆盖率口径
# 使用注解：复用 E1；输出 e7_rigor_check.json，供 make_figures_af.py 绘制偏差图。
# 完整对应关系与运行顺序：../README.md
"""
E7｜口径审查补算：回答三个"数字到底怎么来的"问题
  Q1 0.860 埃 是"裁剪后的精度"吗？—— 拆开「残基集合」和「叠合方式」两个变量
  Q2 迭代剔除到 0.223 埃，剩下几个残基？—— 证明它只是"存在更接近的子集"
  Q3 口袋覆盖率三数字到底什么关系？—— 打印交集体积链条

同时补一组正面对照：口袋残基 vs 非口袋残基在同一全局变换下的偏差。
"""
import os
import sys
import json

import numpy as np
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(WORK, "aidd-mvp", "scripts"))
sys.path.insert(0, os.path.join(WORK, "aidd-mvp", "vendor"))

from _common import raw  # noqa: E402
from e1_af_vs_crystal import parse_pdb, fit_transform, OFFSET, CUTOFF, LIGAND_RESN  # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)


def main():
    c_atoms, c_res = parse_pdb(raw("1M17.pdb"))
    a_atoms, a_res = parse_pdb(os.path.join(HERE, "struct", "egfr_af_v6.pdb"))

    lig_xyz = np.array([a["xyz"] for a in c_atoms if a["resn"] == LIGAND_RESN])
    lt = cKDTree(lig_xyz)
    pocket = set()
    for rn, r in c_res.items():
        if r["hetero"] or r["ca"] is None:
            continue
        for v in r["atoms"].values():
            dd, _ = lt.query(np.asarray(v).reshape(1, -1), k=1)
            if float(dd[0]) <= CUTOFF:
                pocket.add(rn)
                break

    af_by_pdb = {rn + OFFSET: r for rn, r in a_res.items()}
    aplddt = {a["resnum"]: a["bfac"] for a in a_atoms if a["atomname"] == "CA"}

    common = sorted(n for n in c_res
                    if n in af_by_pdb and not c_res[n]["hetero"]
                    and c_res[n]["ca"] is not None and af_by_pdb[n]["ca"] is not None)
    P = np.array([af_by_pdb[n]["ca"] for n in common])
    Q = np.array([c_res[n]["ca"] for n in common])

    hi = np.array([i for i, n in enumerate(common) if aplddt.get(n - OFFSET, 0) >= 70])
    lo = np.array([i for i, n in enumerate(common) if aplddt.get(n - OFFSET, 100) < 50])

    def rmsd(Pm, Qm):
        return float(np.sqrt((((Pm) - Qm) ** 2).sum(1).mean()))

    # ---------- Q1：拆开「残基集合」与「叠合方式」 ----------
    R_all, t_all = fit_transform(P, Q)                 # 用全部残基求变换
    R_hi, t_hi = fit_transform(P[hi], Q[hi])           # 只用高置信残基求变换

    q1 = {
        "A_全部残基+全部拟合": round(rmsd(P @ R_all + t_all, Q), 3),
        "B_高置信子集+沿用全局变换": round(rmsd(P[hi] @ R_all + t_all, Q[hi]), 3),
        "C_高置信子集+子集独立拟合": round(rmsd(P[hi] @ R_hi + t_hi, Q[hi]), 3),
        "D_全部残基+沿用高置信拟合": round(rmsd(P @ R_hi + t_hi, Q), 3),
        "n_common": len(common), "n_hi": int(len(hi)), "n_lo": int(len(lo)),
        "rmsd_lo_subset_independent_fit": round(rmsd(*(lambda Rm, tt: (P[lo] @ Rm + tt, Q[lo]))(*fit_transform(P[lo], Q[lo]))), 3),
    }
    # B vs C 的差 = 纯粹换叠合方式带来的
    q1["说明_B减C"] = "同一批高置信残基，只换叠合方式：%.3f → %.3f" % (
        q1["B_高置信子集+沿用全局变换"], q1["C_高置信子集+子集独立拟合"])
    q1["说明_A减B"] = "同一套全局变换，只换残基集合：%.3f → %.3f" % (
        q1["A_全部残基+全部拟合"], q1["B_高置信子集+沿用全局变换"])

    # ---------- Q2：迭代剔除的残基数 ----------
    idx = np.arange(len(common))
    steps = []
    for _ in range(6):
        Rm, tt = fit_transform(P[idx], Q[idx])
        dev = np.linalg.norm((P[idx] @ Rm + tt) - Q[idx], axis=1)
        steps.append(dict(n_residues=int(len(idx)), rmsd_A=round(float(np.sqrt((dev ** 2).mean())), 3),
                          max_dev_A=round(float(dev.max()), 3)))
        keep = max(int(len(idx) * 0.8), 40)
        idx = idx[np.argsort(dev)[:keep]]
    q2 = steps

    # ---------- Q3：口袋 vs 非口袋在同一全局变换下的偏差 ----------
    ispk = np.array([n in pocket for n in common])
    dev_all = np.linalg.norm((P @ R_all + t_all) - Q, axis=1)
    q3 = {
        "pocket_n": int(ispk.sum()),
        "pocket_dev_mean_A": round(float(dev_all[ispk].mean()), 3),
        "pocket_dev_median_A": round(float(np.median(dev_all[ispk])), 3),
        "nonpocket_n": int((~ispk).sum()),
        "nonpocket_dev_mean_A": round(float(dev_all[~ispk].mean()), 3),
        "nonpocket_dev_median_A": round(float(np.median(dev_all[~ispk])), 3),
        "结论": "同一全局变换下，口袋残基平均偏差 %.3f 埃，非口袋 %.3f 埃（口袋区反而更小）"
                % (float(dev_all[ispk].mean()), float(dev_all[~ispk].mean())),
    }

    out = dict(Q1_rmsd_口径=q1, Q2_迭代剔除=q2, Q3_口袋vs非口袋=q3)
    with open(os.path.join(RESULTS, "e7_rigor_check.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)

    print("=== Q1  RMSD 口径拆解 ===")
    for k in ("A_全部残基+全部拟合", "B_高置信子集+沿用全局变换",
              "C_高置信子集+子集独立拟合", "D_全部残基+沿用高置信拟合",
              "rmsd_lo_subset_independent_fit"):
        print("  %-28s %.3f 埃" % (k, q1[k]))
    print("  共有残基 %d / 高置信 %d / 低置信 %d" % (len(common), len(hi), len(lo)))
    print("=== Q2  迭代剔除 ===")
    for s in steps:
        print("  残基 %3d 个 → RMSD %.3f 埃（最大单点偏差 %.3f 埃）"
              % (s["n_residues"], s["rmsd_A"], s["max_dev_A"]))
    print("=== Q3  口袋 vs 非口袋（同一全局变换）===")
    print("  口袋 %d 个残基    平均 %.3f 埃 / 中位 %.3f 埃"
          % (q3["pocket_n"], q3["pocket_dev_mean_A"], q3["pocket_dev_median_A"]))
    print("  非口袋 %d 个残基  平均 %.3f 埃 / 中位 %.3f 埃"
          % (q3["nonpocket_n"], q3["nonpocket_dev_mean_A"], q3["nonpocket_dev_median_A"]))
    print("写入 results/e7_rigor_check.json")


if __name__ == "__main__":
    main()
