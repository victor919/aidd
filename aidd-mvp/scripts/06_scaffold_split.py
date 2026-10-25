"""06 · 骨架切分 + ECFP 特征矩阵

输入：data/processed/dataset.csv
输出：data/processed/dataset_with_scaffold.csv、X_ecfp.npy、split_indices.json

要点：随机切分下，测试集骨架在训练集里几乎都见过 —— 模型在"默写"。
      骨架切分把重叠清零，指标会明显变差，但那才是真实水平。
"""
import os
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import GroupShuffleSplit

from _common import cfg, ensure_dirs, processed, log, banner, save_json, exists

RDLogger.DisableLog("rdApp.*")
ensure_dirs()
banner("06 · Split + Featurize：先切分，再谈模型")

src = processed("dataset.csv")
if not exists(src):
    raise SystemExit("请先运行 05_standardize.py")
df = pd.read_csv(src)
log(f"读入 {len(df)} 个分子")

# ---------- 1. Bemis-Murcko 骨架 ----------
def scaffold_of(smiles):
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(smiles=smiles, includeChirality=True)
    except Exception:
        return ""

df["scaffold"] = df["smiles"].map(scaffold_of).fillna("")
log(f"唯一骨架 {df['scaffold'].nunique()} 个（{len(df)} 个分子）")

# ---------- 2. ECFP4 ----------
m = cfg["model"]
gen = rdFingerprintGenerator.GetMorganGenerator(
    radius=int(m["fingerprint_radius"]), fpSize=int(m["fingerprint_size"]))


def featurize(smiles):
    mol = Chem.MolFromSmiles(smiles)
    arr = np.zeros((int(m["fingerprint_size"]),), dtype=np.uint8)
    if mol is None:
        return arr
    DataStructs.ConvertToNumpyArray(gen.GetFingerprint(mol), arr)
    return arr


X = np.stack(df["smiles"].map(featurize).values)
y = df["pIC50_median"].values
np.save(processed("X_ecfp.npy"), X)
log(f"ECFP4 矩阵 X {X.shape}，y 范围 [{y.min():.2f}, {y.max():.2f}]")

# ---------- 3. 两种切分对比 ----------
seed = int(cfg["random_seed"])
test_size = float(m["test_size"])

res = {}
for split_name in ["random", "scaffold"]:
    if split_name == "random":
        rs = np.random.RandomState(seed)
        idx = rs.permutation(len(df))
        n_test = int(len(df) * test_size)
        te, tr = idx[:n_test], idx[n_test:]
    else:
        gss = GroupShuffleSplit(n_splits=1, test_size=test_size, random_state=seed)
        tr, te = next(gss.split(df, groups=df["scaffold"].values))

    tr_s = set(df.iloc[tr]["scaffold"])
    te_s = set(df.iloc[te]["scaffold"])
    overlap = len(tr_s & te_s)
    res[split_name] = {"train": tr.tolist(), "test": te.tolist(),
                       "n_train": int(len(tr)), "n_test": int(len(te)),
                       "n_train_scaffolds": len(tr_s),
                       "n_test_scaffolds": len(te_s),
                       "scaffold_overlap": int(overlap)}
    print(f"\n  [{split_name}] train={len(tr)} test={len(te)}")
    print(f"    训练集骨架 {len(tr_s)}，测试集骨架 {len(te_s)}")
    print(f"    测试集骨架在训练集中出现过的：{overlap} / {len(te_s)}")
    if split_name == "scaffold":
        print(f"    -> 重叠清零，模型面对的是真正没见过的新化学骨架")

save_json(res, processed("split_indices.json"))
log("wrote data/processed/split_indices.json")

df.to_csv(processed("dataset_with_scaffold.csv"), index=False)
log("wrote data/processed/dataset_with_scaffold.csv")

banner(f"06 done · 骨架切分下训练集/测试集骨架重叠 = "
       f"{res['scaffold']['scaffold_overlap']}（随机切分为 "
       f"{res['random']['scaffold_overlap']}）")
