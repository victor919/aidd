#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# 文章关联：AI制药细节篇1_靶点与口袋_公众号文章.html
# 脚本职责：实验 E1 / E2：AF 与晶体偏差、pLDDT 与口袋
# 使用注解：读取 struct/egfr_af_v6.pdb 和正式工程 1M17；输出 results/e1_af_vs_crystal.json。
# 完整对应关系与运行顺序：../README.md
"""
E1｜AlphaFold 预测结构 vs 1M17 晶体结构：能不能拿 AF 模型当对接受体？
E2｜pLDDT 能不能当"口袋探测器"？

对照：1M17（X 射线，EGFR 激酶域，含共晶配体 AQ4=erlotinib）
      AF-P00533-F1（AlphaFold DB）
编号：PDB 1M17 残基号 = UniProt P00533 残基号 - 24（由 03_extract_pocket.py 实测校验）
"""
import os
import sys
import json

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.abspath(os.path.join(HERE, "..", ".."))
MVP = os.path.join(WORK, "aidd-mvp")
sys.path.insert(0, os.path.join(MVP, "scripts"))
sys.path.insert(0, os.path.join(MVP, "vendor"))

from _common import raw  # noqa: E402

RESULTS = os.path.join(HERE, "results")
os.makedirs(RESULTS, exist_ok=True)

OFFSET = -24
CUTOFF = 5.0
LIGAND_RESN = "AQ4"


def parse_pdb(path):
    residues = {}
    atoms = []
    seen_models = 0
    with open(path) as f:
        for line in f:
            if line.startswith("MODEL"):
                seen_models += 1
                if seen_models > 1:
                    break
                continue
            if line.startswith("ENDMDL"):
                break
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
            try:
                xyz = np.array([float(line[30:38]), float(line[38:46]),
                                float(line[46:54])], dtype=float)
            except ValueError:
                continue
            bfac = float(line[60:66]) if line[60:66].strip() else 0.0
            atoms.append(dict(resnum=resnum, resn=resname, atomname=atomname,
                              xyz=xyz, bfac=bfac))
            r = residues.setdefault(resnum, dict(resn=resname, atoms={}, hetero=False))
            if atomname not in r["atoms"]:
                r["atoms"][atomname] = xyz
            if line.startswith("HETATM"):
                r["hetero"] = True
    for r in residues.values():
        r["ca"] = r["atoms"].get("CA")
    return atoms, residues


def fit_transform(P, Q):
    """求把 P 叠到 Q 的刚体变换，返回 (Rm, t)，使 P @ Rm + t ≈ Q。

    注意：P/Q 是「行向量」点云（每行一个原子）。Kabsch 得到的真实旋转是
    R = V @ D @ U^T，但行向量点云要写成 (R @ p)^T = p^T @ R^T，
    故直接返回 Rm = R^T = U @ D @ V^T，应用方式就是 X @ Rm + t。
    """
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    Rm = U @ np.diag([1.0, 1.0, d]) @ Vt       # = R^T
    t = Q.mean(axis=0) - P.mean(axis=0) @ Rm
    return Rm, t


def main():
    pdb_crystal = raw("1M17.pdb")
    pdb_af = os.path.join(HERE, "struct", "egfr_af_v6.pdb")

    c_atoms, c_res = parse_pdb(pdb_crystal)
    lig = [a for a in c_atoms if a["resn"] == LIGAND_RESN]
    if not lig:
        raise SystemExit("1M17 中没找到配体 %s" % LIGAND_RESN)
    lig_xyz = np.array([a["xyz"] for a in lig])
    lig_center = lig_xyz.mean(axis=0)

    # 晶体口袋残基：任一重原子距配体任一重原子 <= 5 埃
    # （注意：不能用"距配体质心"——厄洛替尼是细长分子，质心 5 埃内只剩 3 个残基）
    from scipy.spatial import cKDTree as _KDT
    _lt = _KDT(lig_xyz)
    pocket_resnums = set()
    for rnum, r in c_res.items():
        if r["hetero"] or r["ca"] is None:
            continue
        for xyz in r["atoms"].values():
            dd, _ = _lt.query(np.asarray(xyz).reshape(1, -1), k=1)
            if float(dd[0]) <= CUTOFF:
                pocket_resnums.add(rnum)
                break

    a_atoms, a_res = parse_pdb(pdb_af)
    # AF 残基号 = UniProt 编号；PDB 编号 = UniProt 编号 + OFFSET(OFFSET=-24)
    # 故 AF 残基 rn 对应 PDB 编号 rn + OFFSET
    af_by_uniprot = {rn + OFFSET: r for rn, r in a_res.items()}
    af_plddt = {}
    for a in a_atoms:
        if a["atomname"] == "CA":
            af_plddt[a["resnum"]] = a["bfac"]

    common = [n for n in sorted(c_res)
              if n in af_by_uniprot
              and c_res[n]["ca"] is not None
              and af_by_uniprot[n]["ca"] is not None
              and not c_res[n]["hetero"]]
    P = np.array([af_by_uniprot[n]["ca"] for n in common])
    Q = np.array([c_res[n]["ca"] for n in common])
    R, t = fit_transform(P, Q)
    rmsd_global = float(np.sqrt((((P @ R + t) - Q) ** 2).sum(axis=1).mean()))

    # 口袋残基在全局叠合后的偏差
    pocket_common = [n for n in common if n in pocket_resnums]
    Ppk = np.array([af_by_uniprot[n]["ca"] for n in pocket_common]) @ R + t
    Qpk = np.array([c_res[n]["ca"] for n in pocket_common])
    pocket_dev = []
    for n, p in zip(pocket_common, Ppk):
        q = c_res[n]["ca"]
        pocket_dev.append(dict(pdb_resnum=n, uniprot=n - OFFSET,
                               dev_A=round(float(np.linalg.norm(p - q)), 3)))
    pocket_dev_sorted = sorted(pocket_dev, key=lambda x: -x["dev_A"])
    nonpocket_dev = [float(np.linalg.norm(p - q))
                     for n, p in zip(common, P @ R + t)
                     if n not in pocket_resnums
                     for q in [c_res[n]["ca"]]]

    # 口袋区内部结构是否保持：口袋残基两两距离的自变化（配体结合位点的"凹陷"是否还在）
    # 用口袋残基 CA 之间的内部距离矩阵 RMSD（叠合后）
    def internal_dist_matrix(resnums, coords):
        n = len(resnums)
        M = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                d_ = float(np.linalg.norm(coords[i] - coords[j]))
                M[i, j] = M[j, i] = d_
        return M

    pk_idx = {n: i for i, n in enumerate(pocket_common)}
    all_idx = {n: i for i, n in enumerate(common)}
    pk_pos = [all_idx[n] for n in pocket_common]
    M_cry = internal_dist_matrix(pocket_common, Q)
    M_af = internal_dist_matrix(pocket_common, P @ R + t)
    pocket_internal_rmsd = float(np.sqrt(((M_af - M_cry) ** 2).mean()))

    # ---- AF 自己"认为"的口袋：在 AF 坐标里，配体质心变换过去后 5 埃内的残基
    af_lig_center = lig_center @ R + t
    af_res_ca = {n: r["ca"] for n, r in af_by_uniprot.items() if r["ca"] is not None}
    af_near = set()
    for n, c in af_res_ca.items():
        if np.linalg.norm(c - af_lig_center) <= CUTOFF:
            af_near.add(n)
    inter = af_near & pocket_resnums
    union = af_near | pocket_resnums

    # ---- 按 pLDDT 分层看 RMSD（全长模型含 22.8% 无序区，会拉爆全局 RMSD）
    def rmsd_of(mask):
        Rm_, t_ = fit_transform(P[mask], Q[mask])
        return round(float(np.sqrt((((P[mask] @ Rm_ + t_) - Q[mask]) ** 2)
                                .sum(1).mean())), 3)

    hi_mask = [i for i, n in enumerate(common)
               if af_plddt.get(n - OFFSET, 0) >= 70]
    lo_mask = [i for i, n in enumerate(common)
               if af_plddt.get(n - OFFSET, 100) < 50]
    rmsd_hi = rmsd_of(hi_mask) if len(hi_mask) > 20 else None
    rmsd_lo = rmsd_of(lo_mask) if len(lo_mask) > 20 else None

    idx = np.arange(len(common))
    refined = []
    for _ in range(6):
        Rm_, t_ = fit_transform(P[idx], Q[idx])
        dev = np.linalg.norm((P[idx] @ Rm_ + t_) - Q[idx], axis=1)
        refined.append(round(float(np.sqrt((dev ** 2).mean())), 3))
        keep = max(int(len(idx) * 0.8), 40)
        idx = idx[np.argsort(dev)[:keep]]

    # ---- pLDDT 分布
    def stat(v):
        if not v:
            return None
        return dict(n=len(v), mean=round(float(np.mean(v)), 2),
                    median=round(float(np.median(v)), 2),
                    p10=round(float(np.percentile(v, 10)), 2),
                    p90=round(float(np.percentile(v, 90)), 2),
                    min=round(float(np.min(v)), 2),
                    max=round(float(np.max(v)), 2))

    pocket_u = [n - OFFSET for n in pocket_common]
    pocket_plddt = [af_plddt[u] for u in pocket_u if u in af_plddt]
    nonpocket_u = [n - OFFSET for n in common if n not in pocket_resnums]
    nonpocket_plddt = [af_plddt[u] for u in nonpocket_u if u in af_plddt]
    all_plddt = [v for v in af_plddt.values()]

    # 激酶域（PDB 672-995，对应 UniProt 696-1019）与全长的 pLDDT 对比
    kd_u = [u for u in af_plddt if 672 <= (u + OFFSET) <= 995]
    kd_plddt = [af_plddt[u] for u in kd_u]
    out_plddt = [u for u in af_plddt if not (672 <= (u + OFFSET) <= 995)]
    out_domain_plddt = [af_plddt[u] for u in out_plddt]
    frac_low = round(sum(1 for v in all_plddt if v < 50) / len(all_plddt), 4)
    frac_high = round(sum(1 for v in all_plddt if v > 90) / len(all_plddt), 4)

    ranking = [(u, v, (u + OFFSET) in pocket_resnums)
               for u, v in af_plddt.items()]
    ranking.sort(key=lambda x: -x[1])
    topk = {}
    for k in (20, 40, 60, 100, 200):
        sel = ranking[:k]
        hit = sum(1 for x in sel if x[2])
        topk[str(k)] = dict(hit=hit, precision=round(hit / k, 4))
    n_pocket_total = sum(1 for x in ranking if x[2])
    baseline = round(n_pocket_total / len(ranking), 5)

    # 口袋残基在 pLDDT 排序里的名次（1 = 最高置信）
    pocket_ranks = sorted(i + 1 for i, x in enumerate(ranking) if x[2])
    pocket_percentile = [round(100.0 * r / len(ranking), 2) for r in pocket_ranks]

    out = {
        "meta": {
            "structure_crystal": "1M17 (X-ray, EGFR kinase domain, co-crystal ligand AQ4=erlotinib)",
            "structure_predicted": "AF-P00533-F1 (AlphaFold DB v4 model)",
            "uniprot": "P00533",
            "numbering_offset_pdb_minus_uniprot": OFFSET,
            "pocket_cutoff_A": CUTOFF,
        },
        "e1_alpha_vs_crystal": {
            "n_common_residues": len(common),
            "rmsd_ca_global_A": round(rmsd_global, 3),
            "n_pocket_residues_total": len(pocket_resnums),
        "rmsd_ca_pLDDT_ge70_A": rmsd_hi,
        "n_pLDDT_ge70": len(hi_mask),
        "rmsd_ca_pLDDT_lt50_A": rmsd_lo,
        "n_pLDDT_lt50": len(lo_mask),
        "rmsd_iterative_trim": refined,
            "n_pocket_common": len(pocket_common),
            "pocket_dev_mean_A": round(float(np.mean([d_["dev_A"] for d_ in pocket_dev])), 3),
            "pocket_dev_median_A": round(float(np.median([d_["dev_A"] for d_ in pocket_dev])), 3),
            "pocket_dev_max_A": max(d_["dev_A"] for d_ in pocket_dev),
            "nonpocket_dev_mean_A": round(float(np.mean(nonpocket_dev)), 3),
            "nonpocket_dev_median_A": round(float(np.median(nonpocket_dev)), 3),
            "pocket_internal_distmatrix_rmsd_A": round(pocket_internal_rmsd, 3),
            "af_model_pocket_residues_found": len(af_near),
            "af_vs_crystal_pocket_intersection": len(inter),
            "af_vs_crystal_pocket_jaccard": round(len(inter) / len(union), 4) if union else None,
            "worst_pocket_residues": pocket_dev_sorted[:10],
        },
        "e2_plddt_as_pocket_detector": {
            "n_residues_alpha_model": len(af_plddt),
            "n_residues_crystal": len([n for n in c_res if not c_res[n]["hetero"]]),
            "kinase_domain_plddt": stat(kd_plddt),
            "outside_kinase_domain_plddt": stat(out_domain_plddt),
            "frac_plddt_below_50_whole_model": frac_low,
            "frac_plddt_above_90_whole_model": frac_high,
            "pocket_plddt": stat(pocket_plddt),
            "nonpocket_plddt": stat(nonpocket_plddt),
            "whole_model_plddt": stat(all_plddt),
            "topk_precision": topk,
            "random_baseline_precision": baseline,
            "pocket_residue_pLDDT_ranks": pocket_ranks,
            "pocket_residue_pLDDT_percentile": pocket_percentile,
            "note": "precision = 按 pLDDT 降序取前 K 个残基里属于真实口袋(5埃)的比例",
        },
    }

    with open(os.path.join(RESULTS, "e1_af_vs_crystal.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)

    print("=== E1 AlphaFold vs 1M17 ===")
    for k, v in out["e1_alpha_vs_crystal"].items():
        if k != "worst_pocket_residues":
            print("  %-42s %s" % (k, v))
    print("  最偏的口袋残基:")
    for d_ in pocket_dev_sorted[:8]:
        print("    PDB %d / UniProt %d  %.3f 埃" % (d_["pdb_resnum"], d_["uniprot"], d_["dev_A"]))
    print("=== E2 pLDDT 能否当口袋探测器 ===")
    for k, v in out["e2_plddt_as_pocket_detector"].items():
        print("  %-26s %s" % (k, v))
    print("写入 results/e1_af_vs_crystal.json")


if __name__ == "__main__":
    main()
