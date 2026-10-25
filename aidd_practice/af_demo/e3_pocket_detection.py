#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 文章关联：AI制药细节篇1_靶点与口袋_公众号文章.html
# 脚本职责：实验 E3：几何口袋检测对照
# 使用注解：复用 E1 解析与编号常量；输出 e3_pocket_detection.json；pLDDT 不是口袋标签。
# 完整对应关系与运行顺序：../README.md
"""
E3｜把配体删掉，只给蛋白，几何算法能不能把 ATP 口袋找回来？

方法（LIGSITE 式埋藏度 + 连通聚类，纯 numpy/scipy 自写，无外部依赖）：
  1. 蛋白重原子上建 1 埃网格；距任意重原子 <=1.8 埃判为"被占据"
  2. 每个空格点沿 13 根轴线投射线，正负两侧 12 埃内都撞到蛋白则该轴计 1 分
     → 埋藏度 buriedness（0~13），凹腔高、开阔溶剂区低
  3. 埋藏度 >= 阈值的空格点做 6 连通聚类 → 每簇一个候选口袋
  4. 打分 = 体积 x 平均埋藏度

评估口径（避开"整条裂缝并成一个巨块"导致的排名退化）：
  - 真位点中心 = 晶体口袋残基（配体 5 埃内）在该结构自身的 Cα 质心（与坐标系无关）
  - DCC = 预测口袋簇质心 到 真位点中心 的距离（越小越好）
  - 重叠率：仅在有配体的晶体结构上算，与真实配体占据体积的 recall/precision

四组对照：
  A 1M17 晶体结构（删配体）      —— 黄金标准
  B AF-P00533-F1 全长预测结构    —— 含 22.8% 无序区
  C AF 裁到激酶域（PDB 672-995） —— 实际工程中的正确用法
"""
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
from e1_af_vs_crystal import parse_pdb, OFFSET, CUTOFF, LIGAND_RESN  # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

SPACING = 1.0
OCC_R = 1.8
RAY_LEN = 12
BURY_MIN = 9
MIN_CLUSTER = 50
SOLVENTS = {"HOH", "WAT", "SO4", "PO4", "GOL", "EDO", "PEG", "CL", "NA", "MG", "ZN"}

AXES = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, -1, 0), (1, 0, 1),
        (1, 0, -1), (0, 1, 1), (0, 1, -1), (1, 1, 1), (1, 1, -1),
        (1, -1, 1), (1, -1, -1)]


def build_grid(coords, spacing=SPACING, occ_r=OCC_R, ray_len=RAY_LEN):
    lo = coords.min(axis=0) - 6.0
    hi = coords.max(axis=0) + 6.0
    nx, ny, nz = [int(x) for x in (np.ceil((hi - lo) / spacing) + 1)]
    gx, gy, gz = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    pts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], 1) * spacing + lo
    d, _ = cKDTree(coords).query(pts, k=1)
    occ = (d <= occ_r).reshape(nx, ny, nz)
    steps = int(ray_len / spacing)
    buried = np.zeros(occ.shape, dtype=np.int16)
    for ax in AXES:
        hp = np.zeros(occ.shape, bool)
        hn = np.zeros(occ.shape, bool)
        cp = occ.astype(np.float32)
        cn = occ.astype(np.float32)
        for _ in range(steps):
            cp = ndimage.shift(cp, shift=(-ax[0], -ax[1], -ax[2]), order=0, cval=0.0)
            cn = ndimage.shift(cn, shift=(ax[0], ax[1], ax[2]), order=0, cval=0.0)
            hp |= (cp > 0.5)
            hn |= (cn > 0.5)
        buried += (hp & hn).astype(np.int16)
    return dict(occ=occ, buried=buried, lo=lo, spacing=spacing, shape=(nx, ny, nz))


def detect(grid, bury_min=BURY_MIN, min_cluster=MIN_CLUSTER):
    occ, buried, lo, spacing = grid["occ"], grid["buried"], grid["lo"], grid["spacing"]
    free = ~occ
    mask = free & (buried >= bury_min)
    lab, nlab = ndimage.label(mask, structure=ndimage.generate_binary_structure(3, 1))
    pockets = []
    for i in range(1, nlab + 1):
        sel = (lab == i)
        vol = int(sel.sum()) * (spacing ** 3)
        if vol < min_cluster:
            continue
        idx = np.argwhere(sel)
        pockets.append(dict(
            volume_A3=vol,
            mean_buriedness=round(float(buried[sel].mean()), 2),
            score=round(vol * float(buried[sel].mean()), 1),
            centroid=idx.mean(axis=0) * spacing + lo,
            mask=sel,
        ))
    pockets.sort(key=lambda p: -p["score"])
    return pockets


def run(name, pdb_path, true_resnums_pdb, resnum_filter=None, ligand_xyz=None):
    atoms, residues = parse_pdb(pdb_path)
    heavy = []
    for a in atoms:
        if a["atomname"].startswith("H"):
            continue
        if a["resn"] == LIGAND_RESN or a["resn"] in SOLVENTS:
            continue
        if resnum_filter is not None and a["resnum"] not in resnum_filter:
            continue
        heavy.append(a["xyz"])
    heavy = np.array(heavy)

    grid = build_grid(heavy)
    pockets = detect(grid)

    # 真位点中心：口袋残基在该结构中的 Cα 质心
    xyz = [residues[n]["ca"] for n in true_resnums_pdb
           if n in residues and residues[n]["ca"] is not None]
    if not xyz:
        return dict(name=name, error="找不到真口袋残基")
    true_center = np.array(xyz).mean(axis=0)

    for p in pockets:
        p["dcc_A"] = round(float(np.linalg.norm(p["centroid"] - true_center)), 2)

    # 与真实配体占据体积的重叠（仅晶体结构可算）
    overlap = None
    if ligand_xyz is not None:
        occ, buried, lo, spacing = grid["occ"], grid["buried"], grid["lo"], grid["spacing"]
        nx, ny, nz = grid["shape"]
        gx, gy, gz = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
        pts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], 1) * spacing + lo
        d, _ = cKDTree(ligand_xyz).query(pts, k=1)
        true_mask = (d <= 4.0).reshape(nx, ny, nz) & (~occ)
        n_true = int(true_mask.sum())
        rows = []
        for p in pockets[:5]:
            inter = int((p["mask"] & true_mask).sum())
            rows.append(dict(volume_A3=p["volume_A3"], dcc_A=p["dcc_A"],
                             overlap_points=inter,
                             recall=round(inter / n_true, 4) if n_true else None,
                             precision=round(inter / p["volume_A3"] * (spacing ** 3), 4)))
        overlap = dict(n_true_pocket_points=n_true, top5=rows)

    top = [dict(volume_A3=p["volume_A3"], mean_buriedness=p["mean_buriedness"],
                score=p["score"], dcc_A=p["dcc_A"]) for p in pockets[:6]]
    return dict(
        name=name,
        n_protein_atoms=int(len(heavy)),
        n_pockets_detected=len(pockets),
        true_site_center=[round(float(x), 2) for x in true_center],
        top1_volume_A3=pockets[0]["volume_A3"] if pockets else None,
        top1_dcc_A=pockets[0]["dcc_A"] if pockets else None,
        best_dcc_A=min([p["dcc_A"] for p in pockets]) if pockets else None,
        top6_pockets=top,
        ligand_overlap=overlap,
    )


def main():
    pdb_crystal = raw("1M17.pdb")
    pdb_af = os.path.join(HERE, "struct", "egfr_af_v6.pdb")

    c_atoms, c_res = parse_pdb(pdb_crystal)
    lig = [a for a in c_atoms if a["resn"] == LIGAND_RESN]
    lig_xyz = np.array([a["xyz"] for a in lig])
    lig_center = lig_xyz.mean(axis=0)

    # 口袋残基 = 任一重原子距配体任一重原子 <= 5 埃（不能用质心，细长分子会漏）
    lt = cKDTree(lig_xyz)
    true_resnums = set()
    for rnum, r in c_res.items():
        if r["hetero"] or r["ca"] is None:
            continue
        for v in r["atoms"].values():
            dd, _ = lt.query(np.asarray(v).reshape(1, -1), k=1)
            if float(dd[0]) <= CUTOFF:
                true_resnums.add(rnum)
                break

    _, a_res = parse_pdb(pdb_af)
    # AF 文件用 UniProt 编号；PDB 1M17 的 672-995 对应 UniProt 696-1019（OFFSET=-24）
    # 注意：a_res 的 key 本身就是 UniProt 号，不要再加减 OFFSET
    # （旧版写成 {n - OFFSET ...} 会把区间错移到 UniProt 720-1043，N 端少 24 个残基）
    af_kinase_resnums = {n for n in a_res if 672 <= (n + OFFSET) <= 995}

    out = {
        "meta": {
            "method": "LIGSITE 式埋藏度(13轴) + 6连通聚类，自写 numpy/scipy",
            "spacing_A": SPACING, "occupied_radius_A": OCC_R,
            "ray_length_A": RAY_LEN, "buriedness_threshold": BURY_MIN,
            "min_cluster_A3": MIN_CLUSTER,
            "true_site_def": "1M17 共晶配体 AQ4 5 埃内的残基 Cα 质心（各结构自身坐标系）",
            "n_true_pocket_residues": len(true_resnums),
            "offset_AF_uniprot_minus_pdb": -OFFSET,
        },
        "A_crystal_1M17": run("A 1M17 晶体结构（删除共晶配体）", pdb_crystal,
                              true_resnums, ligand_xyz=lig_xyz),
        "B_alphafold_full": run("B AF-P00533-F1 全长预测结构", pdb_af, true_resnums),
        "C_alphafold_kinase_domain": run("C AF 裁到激酶域(PDB 672-995)", pdb_af,
                                         true_resnums, resnum_filter=af_kinase_resnums),
    }

    with open(os.path.join(RESULTS, "e3_pocket_detection.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)

    for k in ("A_crystal_1M17", "B_alphafold_full", "C_alphafold_kinase_domain"):
        r = out[k]
        print("=== %s ===" % r["name"])
        for kk, vv in r.items():
            if kk in ("name", "top6_pockets", "ligand_overlap"):
                continue
            print("  %-30s %s" % (kk, vv))
        if r.get("ligand_overlap"):
            print("  与真实配体体积重叠(top5):")
            for row in r["ligand_overlap"]["top5"]:
                print("    体积%6d 埃^3  DCC %6.2f 埃  recall %.3f  precision %.3f"
                      % (row["volume_A3"], row["dcc_A"], row["recall"], row["precision"]))
        print("  top6 口袋:")
        for i, p in enumerate(r["top6_pockets"], 1):
            print("    #%d 体积%7d 埃^3  埋藏度%.1f  DCC %6.2f 埃"
                  % (i, p["volume_A3"], p["mean_buriedness"], p["dcc_A"]))
    print("写入 results/e3_pocket_detection.json")


if __name__ == "__main__":
    main()
