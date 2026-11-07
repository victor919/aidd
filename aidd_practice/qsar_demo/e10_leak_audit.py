#!/usr/bin/env python
"""e10 · 泄漏审计：随机划分的 0.5095 里，有多少来自"骨架见过"

核心动作：随机划分下把测试集劈成两半——
  A 组：骨架在训练集出现过（"见过"）
  B 组：骨架没出现过（"没见过"）
同一个模型、同一批预测，只分组看误差。如果 A 组明显好于 B 组，
随机划分的高分就有一部分是在奖励"认骨架"。

另做：精确重复、近重复（Tanimoto > 0.95 / 0.9）、骨架组跨集、测量重复。
输出 results/e10_leak_audit.json
"""
import os
import json
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, r2_score

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PROC = os.path.join(ROOT, "aidd-mvp", "data", "processed")
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)


def log(*a):
    print(*a, flush=True)


df = pd.read_csv(os.path.join(PROC, "dataset_with_scaffold.csv"))
X = np.load(os.path.join(PROC, "X_ecfp.npy"))
splits = json.load(open(os.path.join(PROC, "split_indices.json")))
y = df["pIC50_median"].values

out = {}

# ---------- 1. 精确重复 ----------
n = len(df)
uniq = df["smiles"].nunique()
out["exact_duplicate"] = {"rows": int(n), "unique_smiles": int(uniq), "dup_rows": int(n - uniq)}

# ---------- 2. 骨架组跨集 ----------
tr = np.array(splits["random"]["train"])
te = np.array(splits["random"]["test"])
sc = df["scaffold"].values
sc_tr, sc_te = set(sc[tr]), sc[te]
seen_mask = np.array([s in sc_tr for s in sc[te]])
out["random_scaffold"] = {
    "n_test": int(len(te)),
    "seen": int(seen_mask.sum()),
    "unseen": int((~seen_mask).sum()),
    "seen_frac": round(float(seen_mask.mean()), 4),
    "scaffold_groups_in_train": len(sc_tr),
    "scaffold_groups_split_across": int(sum(1 for s in set(sc_te) if s in sc_tr)),
}

# ---------- 3. 同一模型，按"见过/没见过"分组看误差 ----------
rf = RandomForestRegressor(n_estimators=500, random_state=42, n_jobs=-1)
rf.fit(X[tr], y[tr])
p = rf.predict(X[te])
yt = y[te]
err = np.abs(yt - p)
out["decompose_random"] = {
    "all_MAE": round(float(err.mean()), 4),
    "all_R2": round(float(r2_score(yt, p)), 4),
    "seen_MAE": round(float(err[seen_mask].mean()), 4),
    "seen_R2": round(float(r2_score(yt[seen_mask], p[seen_mask])), 4),
    "unseen_MAE": round(float(err[~seen_mask].mean()), 4),
    "unseen_R2": round(float(r2_score(yt[~seen_mask], p[~seen_mask])), 4),
}
log("随机划分分解：", json.dumps(out["decompose_random"], ensure_ascii=False))

# 对照组：骨架划分下同样分组（应当几乎全部 unseen）
tr2 = np.array(splits["scaffold"]["train"])
te2 = np.array(splits["scaffold"]["test"])
sc_tr2 = set(sc[tr2])
seen2 = np.array([s in sc_tr2 for s in sc[te2]])
rf2 = RandomForestRegressor(n_estimators=500, random_state=42, n_jobs=-1)
rf2.fit(X[tr2], y[tr2])
p2 = rf2.predict(X[te2])
err2 = np.abs(y[te2] - p2)
out["scaffold_seen_check"] = {"seen": int(seen2.sum()), "n": int(len(te2)),
                              "MAE": round(float(err2.mean()), 4)}

# ---------- 4. 近重复（只有骨架划分需要查，随机划分必然有） ----------
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fps = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in df["smiles"]]
tr_fps = [fps[i] for i in tr]
mx = np.array([max(DataStructs.BulkTanimotoSimilarity(fps[i], tr_fps)) for i in te])
out["near_dup_random"] = {
    "gt_0.95": int((mx > 0.95).sum()), "gt_0.9": int((mx > 0.9).sum()),
    "gt_0.8": int((mx > 0.8).sum()), "gt_0.7": int((mx > 0.7).sum()),
    "median": round(float(np.median(mx)), 4),
}
# 这些近重复分子的误差 vs 其余
tight = mx > 0.9
if tight.sum() > 0 and (~tight).sum() > 0:
    out["near_dup_random"]["tight_MAE"] = round(float(err[tight].mean()), 4)
    out["near_dup_random"]["rest_MAE"] = round(float(err[~tight].mean()), 4)

# ---------- 5. 测量重复（标签自噪声的规模） ----------
multi = df[df["n_measurements"] > 1]
out["repeats"] = {
    "n_multi": int(len(multi)),
    "median_std": round(float(np.median(multi["pIC50_std"])), 4),
    "mean_std": round(float(multi["pIC50_std"].mean()), 4),
    "q25_std": round(float(np.quantile(multi["pIC50_std"], .25)), 4),
    "q75_std": round(float(np.quantile(multi["pIC50_std"], .75)), 4),
    "frac_std_gt_0.5": round(float((multi["pIC50_std"] > 0.5).mean()), 4),
    "frac_std_gt_1.0": round(float((multi["pIC50_std"] > 1.0).mean()), 4),
    "max_measurements": int(df["n_measurements"].max()),
}

json.dump(out, open(os.path.join(RES, "e10_leak_audit.json"), "w"),
          ensure_ascii=False, indent=2)
log("wrote results/e10_leak_audit.json")
print(json.dumps(out, ensure_ascii=False, indent=2))
