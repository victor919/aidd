"""05 · 数据清洗：标准化 → 去重 → pIC50 → 聚合

输入：data/raw/raw_activity.csv
输出：data/interim/standardized.csv、data/processed/dataset.csv

要点：
  1) MolStandardize 参数名是 skipStandardization，不是 skipStandardize（写错会全军覆没）
  2) isomericSmiles=True —— 不要顺手把手性洗掉
  3) 同一分子多个测量值时用中位数（抗离群），但必须保留测量次数与标准差
"""
import os
import numpy as np
import pandas as pd

from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize

from _common import cfg, ensure_dirs, raw, interim, processed, log, banner, save_json, exists

RDLogger.DisableLog("rdApp.*")
ensure_dirs()
banner("05 · Curation：标准化 → pIC50 → 聚合（比 GNN 更重要的一步）")

src = raw("raw_activity.csv")
if not exists(src):
    raise SystemExit("请先运行 04_fetch_activity.py")
df = pd.read_csv(src)
log(f"读入 {len(df)} 条记录")

# ---------- 1. 保留精确测量值 ----------
a = cfg["activity"]
before = len(df)
df = df[(df["standard_relation"] == a["relation"]) &
        (df["standard_units"] == a["unit"])].copy()
df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
df = df[df["standard_value"].notna() & (df["standard_value"] > 0)]
log(f"只保留 {a['relation']} / {a['unit']}：{before} -> {len(df)}")

df = df[df["canonical_smiles"].notna()]
log(f"去掉无 SMILES：-> {len(df)}")

# ---------- 2. RDKit 标准化 ----------
params = rdMolStandardize.CleanupParameters()
frag_remover = rdMolStandardize.FragmentRemover()          # 去盐/去小碎片


def standardize(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol = rdMolStandardize.Cleanup(mol, params)        # 执行 RDKit Cleanup
        mol = frag_remover.remove(mol)                     # 移除已知碎片；不保证总是只剩一个片段
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, isomericSmiles=True)  # 保留手性
    except Exception:
        return None


df["smiles"] = df["canonical_smiles"].map(standardize)
n_fail = df["smiles"].isna().sum()
df = df[df["smiles"].notna()].copy()
log(f"标准化成功 {len(df)} 条，失败 {n_fail} 条")
if len(df) == 0:
    raise SystemExit("标准化全部失败，请检查 FragmentRemover/参数名")

# ---------- 3. IC50 -> pIC50 ----------
df["pIC50"] = 9 - np.log10(df["standard_value"].astype(float))
print(f"\n  换算示例：IC50 = 100 nM -> pIC50 = 9 - log10(100) = "
      f"{9 - np.log10(100):.1f}")

# ---------- 4. 数据体检：一眼看出有没有脏数据 ----------
print(f"\n  pIC50 分布：min={df['pIC50'].min():.2f}  median={df['pIC50'].median():.2f}  "
      f"max={df['pIC50'].max():.2f}")
weird = df[df["pIC50"] > 11]
print(f"  pIC50 > 11（即 IC50 < 0.01 nM，异常高，需核查）的记录：{len(weird)} 条")
if len(weird):
    top = weird.nlargest(1, "pIC50").iloc[0]
    print(f"    最离谱的一条：pIC50={top['pIC50']:.2f} "
          f"(IC50={10**(9-top['pIC50']):.2e} nM) —— ChEMBL 中待核查的异常记录")
    print(f"    处理方式：MVP 阶段按 pIC50 <= 11 截断（真实项目需回溯 assay 原文）")
n_cut = len(weird)
df = df[df["pIC50"] <= 11].copy()
log(f"截断后 {len(df)} 条")

df[["molecule_chembl_id", "smiles", "pIC50", "assay_chembl_id"]].to_csv(
    interim("standardized.csv"), index=False)
log("wrote data/interim/standardized.csv")

# ---------- 5. 同一分子的多个测量值：中位数聚合 + 保留可信度 ----------
g = df.groupby("smiles")["pIC50"]
agg = pd.DataFrame({
    "pIC50_median": g.median(),
    "n_measurements": g.size(),
    "pIC50_std": g.std(),
    "pIC50_min": g.min(),
    "pIC50_max": g.max(),
}).reset_index()
agg["pIC50_std"] = agg["pIC50_std"].fillna(0.0)
agg = agg.sort_values("pIC50_median", ascending=False).reset_index(drop=True)

multi = (agg["n_measurements"] > 1).sum()
print(f"\n  聚合后唯一分子：{len(agg)}")
print(f"  被重复测量的分子：{multi}（{multi/len(agg):.1%}）")
if multi:
    sub = agg[agg["n_measurements"] > 1]
    print(f"  重复测量的分歧（std 中位数）：{sub['pIC50_std'].median():.2f} 个对数单位")
    print(f"    -> 说明「标签是绝对的」这个假设从数据进来第一天就不成立")

agg.to_csv(processed("dataset.csv"), index=False)
log("wrote data/processed/dataset.csv")

save_json({
    "n_records_in": int(before),
    "n_records_after_filter": int(len(df)),
    "n_standardize_failed": int(n_fail),
    "n_outliers_removed": int(n_cut),
    "n_molecules": int(len(agg)),
    "n_multi_measured": int(multi),
    "median_std_of_multi": float(sub["pIC50_std"].median()) if multi else 0.0,
    "pIC50_range": [float(agg["pIC50_median"].min()), float(agg["pIC50_median"].max())],
}, processed("05_curation.json"))

banner(f"05 done · {len(agg)} 个分子（从 {before} 条记录聚合而来）")
