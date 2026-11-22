# 文章关联：第4课素材（对应正文尚未在相邻目录找到）
# 脚本职责：第4课素材：Vina 与真实活性的对照
# 使用注解：读取正式工程结构和活性数据，执行 Vina；输出 d2 摘要与逐分子 CSV，不运行 DiffDock。
# 完整对应关系与运行顺序：../README.md
"""D2 · Vina 到底能不能排序：用真实活性全面体检

MVP 的 09_docking.py 在 30 个分子上算出了 Vina 分数，但那一轮的结论
是 "pred vs docking Pearson r = -0.059"——**负相关**。

问题：那次只有 30 个分子，且都是按"预测活性高"挑出来的（选择偏差）。
本脚本做一次干净的重测：

问题 1｜在更合理的样本上，Vina 分数与真实 pIC50 相关吗？
       → 从 10468 个有真实活性的分子中随机抽样，跑 Vina，算相关。

问题 2｜Vina 分数能不能区分"活性高"和"活性低"两组？
       → 分箱统计：按 Vina 分数分 5 箱，看每箱的平均真实 pIC50。

问题 3｜如果只能选 10 个分子做实验，按 Vina 分数选 vs 按真实活性选，差多少？
       → 模拟"预算有限"的真实决策场景。

问题 4｜把"预测分数"和"对接分数"混在一起用（consensus），是加分还是减分？
       → 对比三种排序策略在同一预算下的召回。

输出：results/d2_vina_real_activity.json
"""
import os
import sys
import json
import random
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
from meeko.polymer import Polymer, ResidueChemTemplates  # noqa: E402
from meeko import MoleculePreparation  # noqa: E402
from meeko.writer import PDBQTWriterLegacy  # noqa: E402
from vina import Vina  # noqa: E402

SEED = int(cfg["random_seed"])
STRUCT = cfg["structure"]
D = cfg["docking"]
BOX_CENTER = [22.014, 0.253, 52.794]   # 来自 09_docking.json（AQ4 口袋中心）
BOX_SIZE = [float(x) for x in D["box_size"]]
EXHAUST = int(D["exhaustiveness"])

AA20 = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER "
           "THR TRP TYR VAL".split())


def build_receptor():
    """构建受体（去共晶配体 AQ4 与水）——与 09_docking.py 同一套写法"""
    receptor_pdbqt = interim("receptor.pdbqt")
    if os.path.exists(receptor_pdbqt) and os.path.getsize(receptor_pdbqt) > 0:
        print(f"  复用已有受体 {os.path.basename(receptor_pdbqt)}")
        return receptor_pdbqt

    templates = ResidueChemTemplates.create_from_defaults()
    pdb_text = open(raw(f"{STRUCT['pdb_id']}.pdb"), encoding="utf-8", errors="ignore").read()
    to_del = set()
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            resn = line[17:20].strip()
            if resn not in AA20:
                to_del.add(f"{line[21] if len(line) > 21 else 'A'}:{int(line[22:26])}")
    polymer = Polymer.from_pdb_string(
        pdb_text, templates,
        mk_prep=MoleculePreparation(),
        residues_to_delete=sorted(to_del),
        default_altloc=str(D["default_altloc"]),
        allow_bad_res=bool(D["allow_bad_res"]))
    rec_txt = PDBQTWriterLegacy.write_string_from_polymer(polymer)
    rec_txt = rec_txt[0] if isinstance(rec_txt, tuple) else rec_txt
    with open(receptor_pdbqt, "w") as f:
        f.write(rec_txt)
    n_atom = len([l for l in rec_txt.splitlines() if l.startswith("ATOM")])
    print(f"  受体就绪：{n_atom} 个原子，含配体={'AQ4' in rec_txt}")
    return receptor_pdbqt


def prep_ligand(smiles, idx):
    """SMILES -> 3D（带显式 H）-> pdbqt
    meeko 0.8 要求分子带显式氢（'RDKit molecule has implicit Hs'），
    由 meeko 自行处理非极性氢合并；不能自己 RemoveHs。
    """
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
    pdbqt_string, ok, err = PDBQTWriterLegacy.write_string(s)
    if not ok:
        return None
    path = interim(os.path.join("d2_ligs", f"lig{idx:04d}.pdbqt"))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(pdbqt_string)
    return path


def dock(lig_pdbqt, rec_pdbqt):
    v = Vina(sf_name="vina", cpu=int(D["cpu"]), seed=SEED, verbosity=0)
    v.set_receptor(rec_pdbqt)
    v.set_ligand_from_file(lig_pdbqt)
    v.compute_vina_maps(center=BOX_CENTER, box_size=BOX_SIZE)
    v.dock(exhaustiveness=EXHAUST, n_poses=9)
    e = v.energies(n_poses=9)
    return float(e[0][0]), e


def main():
    random.seed(SEED)
    np.random.seed(SEED)

    ds = pd.read_csv(processed("dataset_with_scaffold.csv"))
    print(f"全量数据集：{len(ds)} 个分子")

    # 抽样策略：从有真实活性的分子里分层抽样，保证活性分布覆盖全区间
    ds = ds.sort_values("pIC50_median").reset_index(drop=True)
    n_bins = 8
    ds["_bin"] = pd.cut(ds["pIC50_median"], bins=n_bins, labels=False, include_lowest=True)
    per_bin = 12          # 每箱 12 个 → 96 个分子
    parts = []
    for b, g in ds.groupby("_bin"):
        parts.append(g.sample(min(per_bin, len(g)), random_state=SEED))
    sample = pd.concat(parts).drop(columns=["_bin"]).reset_index(drop=True)
    print(f"分层抽样：{len(sample)} 个分子，活性区间 "
          f"{sample.pIC50_median.min():.2f}–{sample.pIC50_median.max():.2f}")

    rec = build_receptor()
    print(f"受体就绪：{os.path.basename(rec)}")

    rows = []
    err_count = 0
    for i, r in sample.iterrows():
        p = prep_ligand(r["smiles"], i)
        if p is None:
            rows.append({**r.to_dict(), "vina_score": None, "vina_best": None,
                         "error": "prep_failed"})
            err_count += 1
            continue
        try:
            best, energies = dock(p, rec)
            allsc = [float(x[0]) for x in energies]
            rows.append({**r.to_dict(), "vina_score": best, "vina_best": best,
                         "vina_all": allsc, "n_poses": len(allsc)})
        except Exception as e:
            rows.append({**r.to_dict(), "vina_score": None, "vina_best": None,
                         "error": f"{type(e).__name__}: {e}"})
            err_count += 1
            if err_count <= 3:
                print(f"  [失败样例] {r['smiles'][:60]} → {type(e).__name__}: {e}")
        if (i + 1) % 12 == 0:
            ok = sum(1 for x in rows if x.get("vina_score") is not None)
            print(f"  已对接 {i+1}/{len(sample)}（成功 {ok}）")

    out = pd.DataFrame(rows)
    ok = out[out["vina_score"].notna()].copy()
    print(f"\n成功对接 {len(ok)}/{len(out)}")

    res = {"n_total": len(out), "n_docked": len(ok),
           "sampling": f"按 pIC50_median 分 {n_bins} 箱，每箱抽 {per_bin} 个",
           "box_center": BOX_CENTER, "box_size": BOX_SIZE,
           "exhaustiveness": EXHAUST, "seed": SEED}

    if len(ok) >= 10:
        from scipy.stats import spearmanr, pearsonr
        v = ok["vina_score"].values
        y = ok["pIC50_median"].values
        rho_s, p_s = spearmanr(v, y)
        r_p, p_p = pearsonr(v, y)
        res["corr_vina_vs_true_pIC50"] = {
            "spearman_rho": round(float(rho_s), 4), "spearman_p": float(f"{p_s:.3g}"),
            "pearson_r": round(float(r_p), 4), "pearson_p": float(f"{p_p:.3g}"),
            "note": "Vina 分数越负=结合越强，预期与 pIC50 正相关",
        }
        print(f"[问题1] Vina vs 真实 pIC50：Spearman ρ={rho_s:.4f} (p={p_s:.3g})  "
              f"Pearson r={r_p:.4f} (p={p_p:.3g})")

        # 问题 2：分箱
        ok["_vbin"] = pd.qcut(ok["vina_score"], 5, labels=False, duplicates="drop")
        bins = []
        for b, g in ok.groupby("_vbin"):
            bins.append({
                "bin": int(b),
                "vina_range": [round(float(g["vina_score"].min()), 3), round(float(g["vina_score"].max()), 3)],
                "n": len(g),
                "mean_true_pIC50": round(float(g["pIC50_median"].mean()), 3),
                "median_true_pIC50": round(float(g["pIC50_median"].median()), 3),
            })
        res["vina_bins"] = bins
        print("[问题2] 按 Vina 分数分 5 箱后的真实活性：")
        for x in bins:
            print(f"   箱{x['bin']} Vina[{x['vina_range'][0]:.2f},{x['vina_range'][1]:.2f}] "
                  f"n={x['n']:2d}  平均 pIC50={x['mean_true_pIC50']:.3f}")

        # 问题 3/4：有限预算下的召回
        budgets = [10, 20, 30]
        top_true = set(ok.nlargest(30, "pIC50_median")["smiles"])
        strategies = {}
        # 策略 A：只按 Vina 挑（最负的 30 个）
        a30 = set(ok.nsmallest(30, "vina_score")["smiles"])
        # 策略 B：只按真实活性挑
        b30 = set(ok.nlargest(30, "pIC50_median")["smiles"])
        # 策略 C：consensus（活性预测 + 对接）—— 这里用真实活性替代预测活性做上界估计
        ok["_cons"] = 0.7 * ok["pIC50_median"] - 0.3 * ok["vina_score"]
        c30 = set(ok.nlargest(30, "_cons")["smiles"])
        for sname, s30 in [("A_只按Vina", a30), ("B_只按真实活性", b30), ("C_consensus", c30)]:
            rec_map = {}
            for bud in budgets:
                hit = len(set(list(s30)[:bud]) & top_true)
                rec_map[str(bud)] = hit
            strategies[sname] = rec_map
        res["budget_recall_of_top30"] = {
            "top_true_set_size": 30,
            "strategies": strategies,
            "note": "命中数 = 选出的前 N 个里有多少落在真实活性 Top-30 中（上限 30）",
        }
        print("[问题3/4] 有限预算下命中真实活性 Top-30 的个数：")
        for sname, mp in strategies.items():
            print(f"   {sname:16s} " + "  ".join(f"预算{k}→{v2}" for k, v2 in mp.items()))

    ok.to_csv(os.path.join(OUT, "d2_vina_vs_activity.csv"), index=False)
    with open(os.path.join(OUT, "d2_vina_real_activity.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print(f"\n已写出 {OUT}/d2_vina_real_activity.json")


if __name__ == "__main__":
    main()
