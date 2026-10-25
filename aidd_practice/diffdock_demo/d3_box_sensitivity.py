# 文章关联：第4课素材（对应正文尚未在相邻目录找到）
# 脚本职责：第4课素材：正确、偏移与盲搜盒子对照
# 使用注解：依赖 D2 的 CSV；当前目录未发现 d3_box_sensitivity.json，绘图前需完成 D3。
# 完整对应关系与运行顺序：../README.md
"""D3 · 盒子瞄哪，答案就是哪：Vina 分数对空间定义的依赖

D2 已经证明：盒子对准真口袋时，Vina 分数与真实 pIC50 中等相关（rho=-0.4984）。
但那个结论有一个没说出口的前提——你事先知道口袋在哪，而且盒子没瞄错。

本脚本问一个更狠的问题：

    如果不告诉 Vina 口袋在哪，它还能不能给出有用的排序？

做法：同一批有真实活性的分子，在三个盒子里各对接一次。

    盒子 A｜正确口袋  center = AQ4 质心, size 20x20x20
    盒子 B｜偏移 8 埃  沿"口袋->蛋白外"方向平移（模拟大概知道在这附近但瞄偏了）
    盒子 C｜全蛋白盲搜 center = 蛋白包围盒中心, size 覆盖整个蛋白

每个盒子回答两个问题：
    1｜排序还有效吗？→ Spearman rho(Vina 分数, 真实 pIC50)
    2｜分子跑到哪去了？→ 最优 pose 到真实口袋中心的距离

如果盒子 C 的 rho 塌到 0 附近、且 pose 落在离口袋很远的地方，
那就说明：Vina 的分数是"搜索区域的函数"，不是"分子的函数"。
这正是 DiffDock 想解决的那类问题（盲对接、不需要预先定义口袋）。

输出：results/d3_box_sensitivity.json、results/d3_box_detail.csv
"""
import os
import sys
import json
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
MVP = os.path.join(HERE, "..", "..", "aidd-mvp")
sys.path.insert(0, os.path.join(MVP, "scripts"))
sys.path.insert(0, os.path.join(MVP, "vendor"))

OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)

from _common import cfg, raw, interim, processed  # noqa: E402
from meeko import MoleculePreparation  # noqa: E402
from meeko.writer import PDBQTWriterLegacy  # noqa: E402
from vina import Vina  # noqa: E402

SEED = int(cfg["random_seed"])
STRUCT = cfg["structure"]
D = cfg["docking"]
BOX_SIZE = [float(x) for x in D["box_size"]]
N_POSES = int(D["n_poses"])

# 真口袋中心（来自 09_docking.json：共晶配体 AQ4 的 29 个原子的质心）
POCKET_CENTER = np.array([22.014, 0.253, 52.794])

N_MOLECULES = 16          # 三盒对比的样本量（同一个分子集合，三个盒子）
EXHAUST_CMP = 2           # 三盒统一用 2，保证可比（盲搜盒子大，exhaustiveness 相同才公平）

RESULTS = []


# ---------------------------------------------------------------- 受体
def receptor_path():
    p = interim("receptor.pdbqt")
    if not os.path.exists(p):
        raise SystemExit(f"找不到受体 {p}，请先跑 aidd-mvp 的 09_docking.py")
    return p


def receptor_coords(pdbqt):
    xyz = []
    for l in open(pdbqt):
        if l.startswith(("ATOM", "HETATM")):
            xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return np.array(xyz)


# ---------------------------------------------------------------- 配体
def prep_ligand(smiles, tag):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = SEED
    if AllChem.EmbedMolecule(mol, ps) != 0:
        ps.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, ps) != 0:
            return None
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=300)
    except Exception:
        pass
    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    if not setups:
        return None
    s = setups[0]
    s.write_pdbqt = True
    try:
        pdbqt_string, ok, err = PDBQTWriterLegacy.write_string(s)
        if not ok:
            return None
    except Exception:
        out = PDBQTWriterLegacy.write_string(s)
        pdbqt_string = out[0] if isinstance(out, tuple) else out
    path = interim(os.path.join("d3_ligs", f"{tag}.pdbqt"))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(pdbqt_string)
    return path


def dock(lig, rec, center, size, exhaust):
    v = Vina(sf_name="vina", cpu=int(D["cpu"]), seed=SEED, verbosity=0)
    v.set_receptor(rec)
    v.set_ligand_from_file(lig)
    v.compute_vina_maps(center=[float(x) for x in center], box_size=[float(x) for x in size])
    v.dock(exhaustiveness=exhaust, n_poses=N_POSES)
    e = v.energies(n_poses=N_POSES)
    scores = [float(x[0]) for x in e]
    # 取最优 pose 的坐标
    pose_path = interim(os.path.join("d3_ligs", "pose_tmp.pdbqt"))
    v.write_poses(pose_path, n_poses=1, overwrite=True)
    xyz = []
    for l in open(pose_path):
        if l.startswith(("ATOM", "HETATM")):
            xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    xyz = np.array(xyz) if xyz else None
    return scores[0], scores, xyz


def main():
    rec = receptor_path()
    rec_xyz = receptor_coords(rec)
    print(f"受体原子数：{len(rec_xyz)}")

    # ---- 三个盒子 ----
    boxA = {"name": "A_正确口袋", "center": POCKET_CENTER.copy(), "size": BOX_SIZE[:]}

    # B 偏移 8 埃：沿"口袋质心 -> 蛋白质心"的反方向（朝蛋白外侧/溶剂侧）
    prot_center = rec_xyz.mean(axis=0)
    outward = POCKET_CENTER - prot_center
    outward = outward / np.linalg.norm(outward)
    boxB = {"name": "B_偏移8埃", "center": POCKET_CENTER + outward * 8.0, "size": BOX_SIZE[:]}

    # C 全蛋白盲搜：盒子覆盖整个蛋白 + 余量；Vina 网格是立方代价，尺寸设上限
    lo, hi = rec_xyz.min(axis=0), rec_xyz.max(axis=0)
    c_center = (lo + hi) / 2.0
    c_size = [float(min(56.0, round(float(hi[i] - lo[i]) + 6.0, 1))) for i in range(3)]
    boxC = {"name": "C_全蛋白盲搜", "center": c_center, "size": c_size}

    boxes = [boxA, boxB, boxC]
    for b in boxes:
        print(f"  盒子 {b['name']}: center={np.round(b['center'],2).tolist()} size={b['size']}")
    print(f"  蛋白质心={np.round(prot_center,2).tolist()}  "
          f"蛋白包围盒={np.round(lo,1).tolist()}-{np.round(hi,1).tolist()}")

    # ---- 分子集合：从 D2 的 92 个分子里分层抽 24 个（三盒同集）----
    d2csv = os.path.join(OUT, "d2_vina_vs_activity.csv")
    if not os.path.exists(d2csv):
        raise SystemExit(f"找不到 {d2csv}，请先跑 d2_vina_real_activity.py")
    ds = pd.read_csv(d2csv)
    ds = ds.sort_values("pIC50_median").reset_index(drop=True)
    ds["_bin"] = pd.qcut(ds["pIC50_median"], 6, labels=False, duplicates="drop")
    per = int(np.ceil(N_MOLECULES / ds["_bin"].nunique()))
    parts = [g.sample(min(per, len(g)), random_state=SEED) for _, g in ds.groupby("_bin")]
    sample = pd.concat(parts).head(N_MOLECULES).reset_index(drop=True)
    print(f"\n三盒对比样本：{len(sample)} 个分子，"
          f"真实 pIC50 {sample.pIC50_median.min():.2f}-{sample.pIC50_median.max():.2f}")

    for bi, b in enumerate(boxes):
        print(f"\n===== 盒子 {b['name']} =====")
        for i, r in sample.iterrows():
            lig = prep_ligand(r["smiles"], f"{bi}_{i:04d}")
            if lig is None:
                print(f"  [{i+1}/{len(sample)}] 准备失败，跳过")
                continue
            try:
                best, allsc, xyz = dock(lig, rec, b["center"], b["size"], EXHAUST_CMP)
            except Exception as e:
                print(f"  [{i+1}/{len(sample)}] 对接失败：{type(e).__name__}: {e}")
                continue
            dist = None
            if xyz is not None and len(xyz):
                dist = float(np.linalg.norm(xyz.mean(axis=0) - POCKET_CENTER))
            RESULTS.append({
                "box": b["name"], "idx": i, "smiles": r["smiles"],
                "pIC50_median": float(r["pIC50_median"]),
                "vina_best": round(best, 3),
                "n_poses": len(allsc),
                "pose_centroid_dist_to_pocket_A": None if dist is None else round(dist, 2),
            })
            # 每跑完一个盒子就增量落盘，避免长跑中途失败全丢
            df_tmp = pd.DataFrame(RESULTS)
            df_tmp.to_csv(os.path.join(OUT, "d3_box_detail.csv"), index=False)
            if (i + 1) % 8 == 0:
                print(f"  已对接 {i+1}/{len(sample)}  最近分数={best:.3f}  "
                      f"到真口袋距离={dist if dist is None else round(dist,2)} 埃")

    df = pd.DataFrame(RESULTS)
    if df.empty:
        raise SystemExit("没有任何成功对接的结果")

    res = {
        "n_molecules": int(df["idx"].nunique()),
        "exhaustiveness": EXHAUST_CMP,
        "n_poses": N_POSES,
        "seed": SEED,
        "pocket_center": [float(x) for x in POCKET_CENTER],
        "boxes": {b["name"]: {"center": [float(x) for x in b["center"]],
                              "size": b["size"]} for b in boxes},
        "per_box": {},
    }

    from scipy.stats import spearmanr

    for b in boxes:
        sub = df[df["box"] == b["name"]]
        if len(sub) < 8:
            continue
        rho, p = spearmanr(sub["vina_best"].values, sub["pIC50_median"].values)
        d = sub["pose_centroid_dist_to_pocket_A"].dropna()
        res["per_box"][b["name"]] = {
            "n_docked": int(len(sub)),
            "spearman_rho_vina_vs_pIC50": round(float(rho), 4),
            "spearman_p": float(f"{p:.3g}"),
            "mean_best_score": round(float(sub["vina_best"].mean()), 3),
            "best_score_range": [round(float(sub["vina_best"].min()), 3),
                                 round(float(sub["vina_best"].max()), 3)],
            "mean_pose_dist_to_true_pocket": None if d.empty else round(float(d.mean()), 2),
            "min_pose_dist_to_true_pocket": None if d.empty else round(float(d.min()), 2),
            "max_pose_dist_to_true_pocket": None if d.empty else round(float(d.max()), 2),
            "frac_pose_within_6A_of_true_pocket": None if d.empty else
                round(float((d < 6).mean()), 4),
        }
        e = res["per_box"][b["name"]]
        print(f"\n[{b['name']}] n={e['n_docked']}  rho={e['spearman_rho_vina_vs_pIC50']:.4f} "
              f"(p={e['spearman_p']:.3g})  最优分数均值={e['mean_best_score']:.3f}  "
              f"pose 到真口袋距离均值={e['mean_pose_dist_to_true_pocket']} 埃  "
              f"落在真口袋 6 埃内={e['frac_pose_within_6A_of_true_pocket']}")

    # 盒子之间的排名重叠：A 选出的前 8 个，在 B/C 里排第几
    a = df[df["box"] == boxes[0]["name"]].sort_values("vina_best")
    top = list(a["idx"].values[:8])
    cross = {}
    for b in boxes[1:]:
        sub = df[df["box"] == b["name"]].sort_values("vina_best").reset_index(drop=True)
        rank = {int(idx): i + 1 for i, idx in enumerate(sub["idx"].values)}
        rr = [rank.get(int(t), None) for t in top]
        cross[b["name"]] = {
            "top8_from_A": [int(t) for t in top],
            "rank_of_those_in_B": rr,
            "n_still_in_top8": sum(1 for x in rr if x is not None and x <= 8),
            "n_missing": sum(1 for x in rr if x is None),
        }
        print(f"\n盒子 A 的 Top-8 在 {b['name']} 里的排名：{rr}  "
              f"仍在前 8 的有 {cross[b['name']]['n_still_in_top8']} 个")
    res["cross_box_rank_overlap"] = cross

    df.to_csv(os.path.join(OUT, "d3_box_detail.csv"), index=False)
    with open(os.path.join(OUT, "d3_box_sensitivity.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print(f"\n已写出 {OUT}/d3_box_sensitivity.json")


if __name__ == "__main__":
    main()
