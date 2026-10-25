# 文章关联：AI制药实践篇B_实操篇_公众号文章.html
# 脚本职责：早期实验 3：标准化、pIC50、骨架切分及 RF 基线
# 使用注解：读取 step2，输出数据表、指纹及基线摘要；与正式工程的最终指标需区分。
# 完整对应关系与运行顺序：README.md
"""EGFR MVP Step3: data curation -> pIC50 -> scaffold split -> ECFP+RF baseline（本机实跑）"""
import os, json, time, warnings
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import DataStructs
from rdkit.Chem import rdFingerprintGenerator

warnings.filterwarnings("ignore")
RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(BASE, "data", "egfr_mvp")
RAW = os.path.join(DATA, "raw_activity.csv")

print("=" * 72)
print("[1] Load raw ChEMBL activities")
print("=" * 72)
raw = pd.read_csv(RAW)
print("raw rows:", len(raw), "| unique molecules:", raw["molecule_chembl_id"].nunique())
print("standard_relation:", raw["standard_relation"].unique().tolist())
print("standard_units   :", raw["standard_units"].unique().tolist())
print("assay_type       :", raw["assay_type"].unique().tolist())

# 删失/异常值体检（MVP 阶段直接丢弃，但先看看有多少）
print("\n-- value sanity --")
v = pd.to_numeric(raw["standard_value"], errors="coerce")
print("non-numeric:", int(v.isna().sum()))
print("<=0        :", int((v <= 0).sum()))
print("min/max nM :", float(v.min()), float(v.max()))

print("\n" + "=" * 72)
print("[2] Standardization (RDKit MolStandardize)")
print("=" * 72)
_chooser = rdMolStandardize.LargestFragmentChooser()   # 去盐/取最大有机片段

def standardize(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol = rdMolStandardize.Cleanup(mol)          # 去金属、归一化官能团表示
        mol = _chooser.choose(mol)                   # 去掉盐离子、溶剂等小碎片
    except Exception:
        return None
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, isomericSmiles=True)  # 保留手性，不洗掉

t0 = time.time()
uniq = raw[["molecule_chembl_id", "canonical_smiles"]].dropna().drop_duplicates("molecule_chembl_id")
print("unique molecules with smiles:", len(uniq))
uniq = uniq.copy()
uniq["std_smiles"] = uniq["canonical_smiles"].map(standardize)
print(f"standardized in {time.time()-t0:.1f}s; failures:", int(uniq["std_smiles"].isna().sum()))

# 手性是否保留（教学点）
chiral_in = sum(1 for s in uniq["canonical_smiles"] if "@" in s)
chiral_out = sum(1 for s in uniq["std_smiles"].dropna() if "@" in s)
print(f"SMILES containing '@' (chiral spec): before={chiral_in}  after={chiral_out}")

print("\n" + "=" * 72)
print("[3] pIC50 and duplicate-measurement aggregation")
print("=" * 72)
df = raw.dropna(subset=["canonical_smiles"]).copy()
df["value_nM"] = pd.to_numeric(df["standard_value"], errors="coerce")
df = df[(df["value_nM"] > 0)]
df["pIC50"] = 9 - np.log10(df["value_nM"])
print("after value filter:", len(df))

df = df.merge(uniq[["molecule_chembl_id", "std_smiles"]], on="molecule_chembl_id", how="left")
df = df.dropna(subset=["std_smiles"])
print("after standardization:", len(df), "| molecules:", df["std_smiles"].nunique())

g = df.groupby("std_smiles")["pIC50"]
agg = pd.DataFrame({
    "pIC50_median": g.median(),
    "pIC50_mean": g.mean(),
    "pIC50_std": g.std(),
    "n_meas": g.count(),
    "pIC50_min": g.min(),
    "pIC50_max": g.max(),
}).reset_index()
agg["range"] = agg["pIC50_max"] - agg["pIC50_min"]
print(f"\naggregated molecules: {len(agg)}")
print("measurement count distribution:")
print(agg["n_meas"].value_counts().sort_index().head(10).to_string())
print(f"  molecules with >=2 measurements: {(agg['n_meas']>=2).sum()}"
      f" ({100*(agg['n_meas']>=2).mean():.1f}%)")
print(f"  median spread (pIC50_max - pIC50_min) for n>=2: {agg.loc[agg['n_meas']>=2,'range'].median():.2f} log units")
print(f"  molecules with spread > 1.0 log unit: {(agg['range']>1.0).sum()}")

worst = agg.sort_values("range", ascending=False).head(5)
print("\n-- 5 molecules with the largest measurement disagreement --")
print(worst[["std_smiles", "n_meas", "pIC50_median", "pIC50_mean", "pIC50_min", "pIC50_max", "range"]].to_string(index=False))

# mean vs median 差异
d = (agg["pIC50_mean"] - agg["pIC50_median"]).abs()
print(f"\n|mean - median| : median={d.median():.4f}, mean={d.mean():.4f}, max={d.max():.4f}")
print(f"  fraction differing by >0.3 log unit: {(d>0.3).mean()*100:.2f}%")

agg.to_csv(os.path.join(DATA, "dataset.csv"), index=False)
print("\nsaved -> dataset.csv", agg.shape)

print("\n" + "=" * 72)
print("[4] Morgan fingerprint (ECFP4, radius=2, 2048 bits)")
print("=" * 72)
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
def featurize(smiles):
    mol = Chem.MolFromSmiles(smiles)
    arr = np.zeros((2048,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(gen.GetFingerprint(mol), arr)
    return arr

t0 = time.time()
X = np.stack(agg["std_smiles"].map(featurize).values)
y = agg["pIC50_median"].values
print(f"X shape {X.shape}, y shape {y.shape}, featurize {time.time()-t0:.1f}s")
print("pIC50: mean=%.3f std=%.3f min=%.2f max=%.2f" % (y.mean(), y.std(), y.min(), y.max()))
print("bit density: %.4f" % X.mean())

print("\n" + "=" * 72)
print("[5] Scaffold split")
print("=" * 72)
scaffolds = [MurckoScaffold.MurckoScaffoldSmiles(smiles=s, includeChirality=True)
             for s in agg["std_smiles"]]
agg["scaffold"] = scaffolds
print("unique Murcko scaffolds:", agg["scaffold"].nunique(), "from", len(agg), "molecules")
print("singleton scaffolds:", int((agg["scaffold"].value_counts() == 1).sum()))
print("\ntop 5 scaffolds:")
for sc, n in agg["scaffold"].value_counts().head(5).items():
    print(f"  n={n:<5d} {sc}")

from sklearn.model_selection import GroupShuffleSplit
splitter = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
tr_idx, te_idx = next(splitter.split(agg, groups=agg["scaffold"]))
# random split 对照
from sklearn.model_selection import train_test_split
rng = np.random.RandomState(42)
perm = rng.permutation(len(agg))
n_te = len(te_idx)
r_tr, r_te = perm[n_te:], perm[:n_te]

tr_sc = set(agg["scaffold"].iloc[tr_idx]); te_sc = set(agg["scaffold"].iloc[te_idx])
print(f"\nscaffold split: train={len(tr_idx)} test={len(te_idx)}")
print(f"  train scaffolds={len(tr_sc)} test scaffolds={len(te_sc)} overlap={len(tr_sc & te_sc)}")
print(f"random   split: train={len(r_tr)} test={len(r_te)}")
r_tr_sc = set(agg['scaffold'].iloc[r_tr]); r_te_sc = set(agg['scaffold'].iloc[r_te])
print(f"  train scaffolds={len(r_tr_sc)} test scaffolds={len(r_te_sc)} overlap={len(r_tr_sc & r_te_sc)}")

print("\n" + "=" * 72)
print("[6] ECFP4 + RandomForest: random split vs scaffold split")
print("=" * 72)
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, r2_score
from scipy.stats import spearmanr

def run(X, y, tr, te, label, n_estimators=500):
    t0 = time.time()
    m = RandomForestRegressor(n_estimators=n_estimators, random_state=42, n_jobs=-1)
    m.fit(X[tr], y[tr])
    p = m.predict(X[te])
    mae = mean_absolute_error(y[te], p)
    rmse = float(np.sqrt(np.mean((y[te] - p) ** 2)))
    r2 = r2_score(y[te], p)
    sp = spearmanr(y[te], p).statistic
    print(f"  {label:<28s} MAE={mae:.3f}  RMSE={rmse:.3f}  R2={r2:.3f}  Spearman={sp:.3f}  ({time.time()-t0:.1f}s)")
    return dict(label=label, mae=round(mae,3), rmse=round(rmse,3), r2=round(r2,3), spearman=round(sp,3))

res = []
print("  -- median aggregation --")
res.append(run(X, y, r_tr, r_te, "ECFP4+RF / random split"))
res.append(run(X, y, tr_idx, te_idx, "ECFP4+RF / scaffold split"))
y_mean = agg["pIC50_mean"].values
print("  -- mean aggregation (for comparison) --")
res.append(run(X, y_mean, tr_idx, te_idx, "ECFP4+RF / scaffold (mean agg)"))

# 只看单点测量分子 vs 全部
mask_single = (agg["n_meas"] == 1).values
print(f"\n  -- sanity: molecules with a single measurement only (n={mask_single.sum()}) --")

# 简单 baseline：预测训练集均值
print("  -- trivial baseline: predict train median --")
for label, tr, te in [("random", r_tr, r_te), ("scaffold", tr_idx, te_idx)]:
    p = np.full(len(te), y[tr].mean())
    print(f"  {label:<28s} MAE={mean_absolute_error(y[te], p):.3f}  R2={r2_score(y[te], p):.3f}")

json.dump(res, open(os.path.join(DATA, "baseline_results.json"), "w"), indent=2)
agg.to_csv(os.path.join(DATA, "dataset_with_scaffold.csv"), index=False)
np.save(os.path.join(DATA, "X_ecfp.npy"), X)
print("\nsaved -> dataset_with_scaffold.csv / X_ecfp.npy / baseline_results.json")
