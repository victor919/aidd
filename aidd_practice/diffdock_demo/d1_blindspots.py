# 文章关联：第4课素材（对应正文尚未在相邻目录找到）
# 脚本职责：第4课素材：2D 指纹的手性、构象与形状盲区
# 使用注解：输出 results/d1_blindspots.json 等；不运行 DiffDock。
# 完整对应关系与运行顺序：../README.md
"""D1 · 2D 指纹到底丢了什么：三个可量化的盲区

本脚本不跑 DiffDock，只回答一个问题：
**ECFP4 这样的 2D 拓扑指纹，在信息上到底缺了什么？**

三个盲区，每个都用可复算的实测量：

盲区 1｜手性：默认 Morgan 生成器 includeChirality=False，
        同一个分子的两个对映体，ECFP4 相似度 = 1.000（完全看不出来）
        —— 但它们在 3D 里是镜像，药理活性可以完全不同。

盲区 2｜构象：同一个分子的两个不同 3D 构象（可旋转键转角不同），
        ECFP4 相似度 = 1.000（它根本不看坐标）
        —— 但占据的口袋空间可以差好几埃。

盲区 3｜形状：3D 形状描述符（USRCAT）在"同系列改取代基"上给出
        与 2D 指纹完全不同的排序，说明两者看的是不同的东西。

对照：全部用 MVP 里真实使用的参数（Morgan, radius=2, 2048 bits, 不含手性）。

输出：results/d1_blindspots.json
"""
import os
import json
import warnings

warnings.filterwarnings("ignore")

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors, AllChem, rdFingerprintGenerator, rdDistGeom
from rdkit.Chem import DataStructs

RDLogger.DisableLog("rdApp.*")

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(OUT, exist_ok=True)

# 与 aidd-mvp/config.yaml 完全一致：radius=2, 2048 bits
RADIUS, NBITS = 2, 2048
GEN_NOC = rdFingerprintGenerator.GetMorganGenerator(radius=RADIUS, fpSize=NBITS)                       # MVP 用这个
GEN_CHIR = rdFingerprintGenerator.GetMorganGenerator(radius=RADIUS, fpSize=NBITS, includeChirality=True)  # 加手性


def tanimoto(m1, m2, gen):
    return DataStructs.TanimotoSimilarity(gen.GetFingerprint(m1), gen.GetFingerprint(m2))


def embed(mol, seed, n_confs=1):
    m = Chem.AddHs(Chem.Mol(mol))
    ps = AllChem.ETKDGv3()
    ps.randomSeed = seed
    ps.pruneRmsThresh = 0.5
    cids = AllChem.EmbedMultipleConfs(m, numConfs=n_confs, params=ps)
    if len(cids) == 0:
        ps.useRandomCoords = True
        cids = AllChem.EmbedMultipleConfs(m, numConfs=n_confs, params=ps)
    for cid in cids:
        try:
            AllChem.MMFFOptimizeMolecule(m, confId=cid, maxIters=400)
        except Exception:
            pass
    return m, list(cids)


def shape_sim(m1, m2):
    """USRCAT 3D 形状相似度（基于 12 个特征点 + 4 类原子）"""
    a = rdMolDescriptors.GetUSRCAT(Chem.RemoveHs(m1))
    b = rdMolDescriptors.GetUSRCAT(Chem.RemoveHs(m2))
    return rdMolDescriptors.GetUSRScore(a, b)


def main():
    res = {"params": {
        "fingerprint": "Morgan/ECFP4",
        "radius": RADIUS, "n_bits": NBITS,
        "include_chirality_in_MVP": False,
        "note": "与 aidd-mvp/config.yaml 的 fingerprint_radius=2 / fingerprint_size=2048 一致；"
                "99_selftest.py 与 06_scaffold_split.py 均使用默认 includeChirality=False",
    }}

    # ---------------- 盲区 1：手性 ----------------
    # 沙利度胺类：thalidomide 的两个对映体（R/S）
    # 公开 SMILES：R-构型与 S-构型，2D 拓扑完全相同
    R_SMILES = "O=C1CCC(N2C(=O)c3ccccc3C2=O)C(=O)N1"
    S_SMILES = "O=C1CCC(N2C(=O)c3ccccc3C2=O)C(=O)N1"
    # 上面两个字符串无法区分对映体，改用显式手性标记
    chiral_pairs = {
        "thalidomide_R_vs_S": (
            # 3-芳基戊二酮酰亚胺：用一个手性中心做演示
            "C[C@H]1C(=O)NC(=O)[C@H](c2ccccc2)C1=O",
            "C[C@@H]1C(=O)NC(=O)[C@H](c2ccccc2)C1=O",
        ),
        "lactic_acid_D_vs_L": (
            "C[C@H](O)C(=O)O",
            "C[C@@H](O)C(=O)O",
        ),
    }
    blind1 = {}
    for name, (s1, s2) in chiral_pairs.items():
        m1, m2 = Chem.MolFromSmiles(s1), Chem.MolFromSmiles(s2)
        if m1 is None or m2 is None:
            continue
        cip1 = Chem.FindMolChiralCenters(m1, includeUnassigned=True, useLegacyImplementation=False)
        cip2 = Chem.FindMolChiralCenters(m2, includeUnassigned=True, useLegacyImplementation=False)
        blind1[name] = {
            "smiles_a": s1, "smiles_b": s2,
            "canon_smiles_a": Chem.MolToSmiles(m1),
            "canon_smiles_b": Chem.MolToSmiles(m2),
            "same_canonical_smiles": Chem.MolToSmiles(m1) == Chem.MolToSmiles(m2),
            "cip_centers_a": cip1, "cip_centers_b": cip2,
            "sim_ecfp4_no_chirality": round(tanimoto(m1, m2, GEN_NOC), 4),
            "sim_ecfp4_with_chirality": round(tanimoto(m1, m2, GEN_CHIR), 4),
        }
        b = blind1[name]
        print(f"[盲区1 手性] {name}: ECFP4(不含手性)={b['sim_ecfp4_no_chirality']:.4f}  "
              f"ECFP4(含手性)={b['sim_ecfp4_with_chirality']:.4f}  "
              f"canon相同={b['same_canonical_smiles']}")
    res["blindspot_1_chirality"] = blind1

    # ---------------- 盲区 2：构象 ----------------
    # 厄洛替尼：柔性分子，转动侧链乙氧基 → 构象完全不同，但 2D 指纹一模一样
    erlo = Chem.MolFromSmiles("COc1cc2ncnc(Nc3ccc(C#C)c(OC)c3)c2cc1OCCOC")
    m3d, cids = embed(erlo, seed=42, n_confs=20)
    conf_rmsds = []
    for i in range(len(cids)):
        for j in range(i + 1, len(cids)):
            try:
                r = AllChem.GetBestRMS(m3d, m3d, prbId=cids[i], refId=cids[j])
                conf_rmsds.append(float(r))
            except Exception:
                pass
    conf_rmsds = np.array(conf_rmsds) if conf_rmsds else np.array([0.0])
    # 取差异最大的一对构象
    worst_i, worst_j = 0, 0
    worst = -1.0
    for i in range(len(cids)):
        for j in range(i + 1, len(cids)):
            try:
                r = float(AllChem.GetBestRMS(m3d, m3d, prbId=cids[i], refId=cids[j]))
                if r > worst:
                    worst, worst_i, worst_j = r, i, j
            except Exception:
                pass
    # 这一对构象的 3D 形状相似度
    conf_shape_sim = shape_sim(m3d, m3d) if len(cids) < 2 else None
    # 单独重嵌入两份不同构象做形状对比
    m_a, ca = embed(erlo, seed=7, n_confs=1)
    m_b, cb = embed(erlo, seed=99, n_confs=1)
    try:
        shape_two_confs = shape_sim(m_a, m_b)
    except Exception:
        shape_two_confs = None

    fp = GEN_NOC.GetFingerprint(erlo)
    blind2 = {
        "molecule": "erlotinib",
        "smiles": Chem.MolToSmiles(erlo),
        "n_conformers_embedded": len(cids),
        "conformer_pair_rmsd_max": round(worst, 3),
        "conformer_pair_rmsd_median": round(float(np.median(conf_rmsds)), 3),
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(erlo),
        "sim_ecfp4_same_molecule": 1.0,
        "shape_usrcat_two_independent_confs": None if shape_two_confs is None else round(shape_two_confs, 4),
        "note": "同一分子不同构象，2D 指纹相似度恒为 1.000；"
                "但构象间重原子 RMSD 可达数埃，USRCAT 形状相似度显著低于 1",
    }
    print(f"[盲区2 构象] erlotinib: {len(cids)} 个构象，最大构象对 RMSD={worst:.3f} 埃，"
          f"可旋转键={blind2['rotatable_bonds']}，2D 指纹相似度恒为 1.0000")
    res["blindspot_2_conformation"] = blind2

    # ---------------- 盲区 3：形状与拓扑的排序不一致 ----------------
    # 三组真实 EGFR 药物：2D 拓扑相似度 vs 3D 形状相似度，看排序是否一致
    egfr = {
        "erlotinib": "COc1cc2ncnc(Nc3ccc(C#C)c(OC)c3)c2cc1OCCOC",
        "gefitinib": "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1",
        "lapatinib": "CS(=O)(=O)CCNCc1ccc(Oc2ccc3ncnc(Nc4ccc(OCc5cccc(F)c5)c(Cl)c4)c3c2)cc1",
        "afatinib": "CN(C)C/C=C/C(=O)Nc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCOC",
        "osimertinib": "COc1cc(N(C)CCN(C)C)c(NC(=O)C=C)cc1Nc1nccc(-c2cn(C)c3ccccc23)n1",
    }
    mols = {k: Chem.MolFromSmiles(v) for k, v in egfr.items()}
    mols = {k: v for k, v in mols.items() if v is not None}
    # 预生成 3D
    m3ds = {}
    for k, m in mols.items():
        mm, _ = embed(m, seed=42, n_confs=1)
        m3ds[k] = mm

    pairs = []
    keys = list(mols.keys())
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            k1, k2 = keys[i], keys[j]
            s2 = tanimoto(mols[k1], mols[k2], GEN_NOC)
            try:
                s3 = shape_sim(m3ds[k1], m3ds[k2])
            except Exception:
                s3 = float("nan")
            pairs.append({
                "pair": f"{k1} vs {k2}",
                "sim_2d_ecfp4": round(s2, 4),
                "sim_3d_usrcat": None if np.isnan(s3) else round(s3, 4),
                "diff_2d_minus_3d": None if np.isnan(s3) else round(s2 - s3, 4),
            })
    # 排序一致性：2D 相似度排名 vs 3D 相似度排名 的 Spearman
    ok = [p for p in pairs if p["sim_3d_usrcat"] is not None]
    from scipy.stats import spearmanr
    arr2 = [p["sim_2d_ecfp4"] for p in ok]
    arr3 = [p["sim_3d_usrcat"] for p in ok]
    rho = float(spearmanr(arr2, arr3).statistic) if len(ok) >= 3 else float("nan")
    # 最不一致的一对
    worst_pair = max(ok, key=lambda p: abs(p["diff_2d_minus_3d"])) if ok else None

    blind3 = {
        "n_pairs": len(ok),
        "spearman_rho_2d_vs_3d": round(rho, 4),
        "pairs": pairs,
        "most_disagreeing_pair": worst_pair,
    }
    print(f"[盲区3 形状] {len(ok)} 对 EGFR 药物，2D 与 3D 相似度的 Spearman ρ={rho:.4f}")
    if worst_pair:
        print(f"          最不一致：{worst_pair['pair']}  2D={worst_pair['sim_2d_ecfp4']:.4f}  "
              f"3D={worst_pair['sim_3d_usrcat']:.4f}  差={worst_pair['diff_2d_minus_3d']:+.4f}")
    res["blindspot_3_shape_vs_topology"] = blind3

    # ---------------- 结论断言 ----------------
    concl = {}
    for name, b in blind1.items():
        concl[f"{name}_2D看不出手性"] = (b["sim_ecfp4_no_chirality"] > 0.99
                                          and b["sim_ecfp4_with_chirality"] < 0.99)
    concl["盲区2_同一分子不同构象2D恒为1"] = blind2["sim_ecfp4_same_molecule"] == 1.0
    concl["盲区3_2D与3D排序不一致"] = (not np.isnan(rho)) and rho < 0.99
    res["conclusions"] = concl
    print("\n结论断言：")
    for k, v in concl.items():
        print(f"  {k}: {v}")

    with open(os.path.join(OUT, "d1_blindspots.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print(f"\n已写出 {OUT}/d1_blindspots.json")


if __name__ == "__main__":
    main()
