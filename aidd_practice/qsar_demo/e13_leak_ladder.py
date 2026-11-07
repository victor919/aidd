#!/usr/bin/env python
"""e13 · 92 分子闭环：把"近邻泄漏"一层层剥掉

e9 里 QSAR 在 92 个分子上 ρ = 0.895，高得可疑——因为这 92 个分子取自同一批 10468，
只把它们本身从训练集删掉，它们的类似物还在里面。这个脚本做一次"泄漏阶梯"：

  A  只删这 92 个分子                    （e9 的做法）
  B  再删掉与它们同骨架的所有训练分子      （骨架级隔离）
  C  再删掉与它们 Tanimoto >= 0.7 的训练分子（近邻级隔离）

同一批 92 个分子、同一个模型、同一个目标集（真实 pIC50 Top-30），
看命中数与秩相关怎么随隔离强度下降。Vina 那一行来自 ⑩ 的既有结果，不重跑。

输出 results/e13_leak_ladder.json
"""
import os
import json
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error
from scipy.stats import spearmanr

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PROC = os.path.join(ROOT, "aidd-mvp", "data", "processed")
D2 = os.path.join(ROOT, "aidd_practice", "diffdock_demo", "results", "d2_vina_vs_activity.csv")
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)

df = pd.read_csv(os.path.join(PROC, "dataset_with_scaffold.csv"))
X = np.load(os.path.join(PROC, "X_ecfp.npy"))
y = df["pIC50_median"].values
smis = df["smiles"].values
sc = df["scaffold"].values

d2 = pd.read_csv(D2)
test_smis = list(d2["smiles"])
yt = d2["pIC50_median"].values
vina = d2["vina_score"].values

idx_map = {s: i for i, s in enumerate(smis)}
te_idx = np.array([idx_map[s] for s in test_smis])
is_test = np.zeros(len(df), bool)
is_test[te_idx] = True

# 相似度矩阵：92 x 10468
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fps = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in smis]
test_fps = [fps[i] for i in te_idx]
S = np.zeros((len(te_idx), len(df)))
for r, f in enumerate(test_fps):
    S[r] = DataStructs.BulkTanimotoSimilarity(f, fps)
max_to_test = S.max(axis=0)          # 每个训练分子离 92 个测试分子最近有多近

target = set(np.argsort(-yt)[:30])


def evaluate(keep_mask, tag):
    tr = np.where(keep_mask)[0]
    rf = RandomForestRegressor(n_estimators=500, random_state=42, n_jobs=-1)
    rf.fit(X[tr], y[tr])
    p = rf.predict(X[te_idx])
    rho = float(spearmanr(yt, p).statistic)
    order_q = np.argsort(-p)
    order_v = np.argsort(vina)          # Vina 越低越好
    res = {"tag": tag, "n_train": int(len(tr)),
           "removed": int(len(df) - len(tr)),
           "MAE": round(float(mean_absolute_error(yt, p)), 4),
           "spearman": round(rho, 4),
           "recall": {str(b): int(len(set(order_q[:b]) & target)) for b in (10, 20, 30)},
           "recall_vina": {str(b): int(len(set(order_v[:b]) & target)) for b in (10, 20, 30)},
           "test_maxsim_median": round(float(np.median(S.max(axis=1))), 4)}
    print(res, flush=True)
    return res, p


rows = []
# A 只删 92
rows.append(evaluate(~is_test, "A 只删这 92 个")[0])

# B 再删同骨架
test_sc = set(sc[te_idx])
mB = (~is_test) & (~np.isin(sc, list(test_sc)))
rows.append(evaluate(mB, "B 再删同骨架")[0])

# C 再删近邻 >= 0.7
mC = mB & (max_to_test < 0.7)
rows.append(evaluate(mC, "C 再删近邻 0.7")[0])

# 每个测试分子在其训练集里的最近邻（用 C 的训练集）
trC = np.where(mC)[0]
nn = np.array([max(DataStructs.BulkTanimotoSimilarity(fps[i], [fps[j] for j in trC]))
               for i in te_idx])

out = {
    "n_test": len(te_idx),
    "ladder": rows,
    "test_nn_vs_trainC": {
        "median": round(float(np.median(nn)), 4),
        "q25": round(float(np.quantile(nn, .25)), 4),
        "q75": round(float(np.quantile(nn, .75)), 4),
        "frac_gt_0.7": round(float((nn > 0.7).mean()), 4),
        "frac_gt_0.5": round(float((nn > 0.5).mean()), 4),
        "frac_lt_0.4": round(float((nn < 0.4).mean()), 4),
    },
    "vina_spearman_raw": round(float(spearmanr(yt, vina).statistic), 4),
    "note": "目标集 = 92 个分子中真实 pIC50 Top-30；所有排序只用预测值或 Vina 分数，无真值参与",
}
json.dump(out, open(os.path.join(RES, "e13_leak_ladder.json"), "w"),
          ensure_ascii=False, indent=2)
print(json.dumps(out, ensure_ascii=False, indent=2))
