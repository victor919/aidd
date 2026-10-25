"""07 · Baseline：ECFP4 + Random Forest

输入：data/processed/dataset_with_scaffold.csv、X_ecfp.npy、split_indices.json
输出：models/rf_<split>.joblib、data/processed/07_baseline.json、predictions.csv

要点：不要一上来就 GNN。先有一个无聊的 baseline，
      否则你不知道高级模型到底"额外获得了什么"。
"""
import os
import json
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, r2_score

from _common import cfg, P, ensure_dirs, processed, path, log, banner, save_json, exists

ensure_dirs()
banner("07 · Baseline：ECFP4 + Random Forest（两种切分对照）")

df = pd.read_csv(processed("dataset_with_scaffold.csv"))
X = np.load(processed("X_ecfp.npy"))
y = df["pIC50_median"].values
splits = json.load(open(processed("split_indices.json")))
log(f"数据 {X.shape}，标签范围 [{y.min():.2f}, {y.max():.2f}]")

m = cfg["model"]
out = {}
for name in ["random", "scaffold"]:
    tr = np.array(splits[name]["train"]); te = np.array(splits[name]["test"])
    rf = RandomForestRegressor(n_estimators=int(m["rf_n_estimators"]),
                               random_state=int(cfg["random_seed"]), n_jobs=-1)
    rf.fit(X[tr], y[tr])
    pred = rf.predict(X[te])
    mae = float(mean_absolute_error(y[te], pred))
    rmse = float(np.sqrt(((y[te] - pred) ** 2).mean()))
    r2 = float(r2_score(y[te], pred))
    out[name] = {"n_train": int(len(tr)), "n_test": int(len(te)),
                 "MAE": round(mae, 4), "RMSE": round(rmse, 4), "R2": round(r2, 4),
                 "scaffold_overlap": splits[name]["scaffold_overlap"]}
    print(f"\n  [{name}] MAE={mae:.4f}  RMSE={rmse:.4f}  R2={r2:.4f}")
    import joblib
    joblib.dump(rf, path(P["models"], f"rf_{name}.joblib"))

    # 第一版预测结果（后续 09/10 要用）
    sub = df.iloc[te][["smiles"]].copy()
    sub["true_pIC50"] = y[te]
    sub["pred_pIC50"] = pred
    sub["split"] = name
    sub.to_csv(processed(f"predictions_{name}.csv"), index=False)

delta = out["scaffold"]["MAE"] - out["random"]["MAE"]
print(f"\n  随机切分 vs 骨架切分：MAE {out['random']['MAE']:.3f} -> "
      f"{out['scaffold']['MAE']:.3f}（+{delta/out['random']['MAE']:.1%}）")
print(f"  这不是模型变差了，是此前的好成绩里有「默写」的成分。")

save_json(out, processed("07_baseline.json"))
log("wrote data/processed/07_baseline.json")
banner(f"07 done · scaffold split MAE = {out['scaffold']['MAE']}")
