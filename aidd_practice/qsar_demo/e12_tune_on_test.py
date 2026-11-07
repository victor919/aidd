#!/usr/bin/env python
"""e12 · 测试集调参：偷看考卷能白拿多少分

三种选超参的方式，都在骨架划分上做：
  A 默认参数（500 树 / max_features=sqrt / leaf=1）
  B 用验证集选（训练集内部再切 20%，按骨架分组）
  C 用测试集选（反复看测试 R²，挑最好的那个）

同一个模型族、同一份训练数据，只改"用什么数据决定超参"，看测试成绩差多少。
输出 results/e12_tune_on_test.json
"""
import os
import json
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_absolute_error
from sklearn.model_selection import GroupShuffleSplit

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PROC = os.path.join(ROOT, "aidd-mvp", "data", "processed")
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)

df = pd.read_csv(os.path.join(PROC, "dataset_with_scaffold.csv"))
X = np.load(os.path.join(PROC, "X_ecfp.npy"))
splits = json.load(open(os.path.join(PROC, "split_indices.json")))
y = df["pIC50_median"].values
sc = df["scaffold"].values

tr = np.array(splits["scaffold"]["train"])
te = np.array(splits["scaffold"]["test"])

# 训练集内部切验证集（按骨架分组，避免验证集自身泄漏）
gss = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
sub_tr, sub_va = next(gss.split(tr, groups=sc[tr]))
tr_fit = tr[sub_tr]
va = tr[sub_va]

GRID = [
    {"n_estimators": 500, "max_features": "sqrt", "min_samples_leaf": 1},   # 默认
    {"n_estimators": 200, "max_features": "sqrt", "min_samples_leaf": 1},
    {"n_estimators": 1000, "max_features": "sqrt", "min_samples_leaf": 1},
    {"n_estimators": 500, "max_features": 0.3, "min_samples_leaf": 1},
    {"n_estimators": 500, "max_features": 0.5, "min_samples_leaf": 1},
    {"n_estimators": 500, "max_features": "sqrt", "min_samples_leaf": 3},
    {"n_estimators": 500, "max_features": "sqrt", "min_samples_leaf": 5},
    {"n_estimators": 1000, "max_features": 0.5, "min_samples_leaf": 3},
]

rows = []
for cfg in GRID:
    rf = RandomForestRegressor(random_state=42, n_jobs=-1, **cfg)
    rf.fit(X[tr_fit], y[tr_fit])
    pv = rf.predict(X[va])
    pt = rf.predict(X[te])
    rows.append({
        **cfg,
        "val_R2": round(float(r2_score(y[va], pv)), 4),
        "val_MAE": round(float(mean_absolute_error(y[va], pv)), 4),
        "test_R2": round(float(r2_score(y[te], pt)), 4),
        "test_MAE": round(float(mean_absolute_error(y[te], pt)), 4),
    })
    print(rows[-1], flush=True)

by_val = max(rows, key=lambda r: r["val_R2"])
by_test = max(rows, key=lambda r: r["test_R2"])
default = rows[0]

out = {
    "n_train_fit": int(len(tr_fit)), "n_val": int(len(va)), "n_test": int(len(te)),
    "grid": rows,
    "default": default,
    "best_by_val": by_val,
    "best_by_test": by_test,
    "inflation_R2": round(by_test["test_R2"] - by_val["test_R2"], 4),
    "inflation_MAE": round(by_val["test_MAE"] - by_test["test_MAE"], 4),
    "val_test_rank_disagree": bool(by_val["n_estimators"] != by_test["n_estimators"]
                                  or by_val["max_features"] != by_test["max_features"]
                                  or by_val["min_samples_leaf"] != by_test["min_samples_leaf"]),
}
json.dump(out, open(os.path.join(RES, "e12_tune_on_test.json"), "w"),
          ensure_ascii=False, indent=2)
print(json.dumps({k: v for k, v in out.items() if k != "grid"}, ensure_ascii=False, indent=2))
