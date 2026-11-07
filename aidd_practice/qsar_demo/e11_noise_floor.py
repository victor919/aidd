#!/usr/bin/env python
"""e11 · 噪声底：模型误差离"实验自己抖多少"还有多远

思路
  同一个分子被不同 assay 反复测，测量值的离散度就是标签的不可约噪声。
  3710 个分子有多次测量，中位标准差 0.4302 个 log 单位。
  把 0.4302 当作单次测量的噪声 sigma，则：
    - 最好的模型能拿到的 R² 上限 = 1 - sigma^2 / var(y)
    - 最好的模型能拿到的 MAE 下限 ≈ sigma * sqrt(2/pi)
  这两个数，就是"0.5095 到底算好还是算差"的分母。

输出 results/e11_noise_floor.json
"""
import os
import json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PROC = os.path.join(ROOT, "aidd-mvp", "data", "processed")
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)

df = pd.read_csv(os.path.join(PROC, "dataset_with_scaffold.csv"))
y = df["pIC50_median"].values
var_y = float(np.var(y, ddof=1))

multi = df[df["n_measurements"] > 1]
sigma_med = float(np.median(multi["pIC50_std"]))
sigma_mean = float(multi["pIC50_std"].mean())
sigma_q75 = float(np.quantile(multi["pIC50_std"], 0.75))

out = {
    "n": int(len(df)),
    "var_y": round(var_y, 4),
    "sd_y": round(float(np.sqrt(var_y)), 4),
    "n_multi": int(len(multi)),
    "sigma_median": round(sigma_med, 4),
    "sigma_mean": round(sigma_mean, 4),
    "sigma_q75": round(sigma_q75, 4),
    "ceilings": {},
    "observed": {"random_MAE": 0.5095, "scaffold_MAE": 0.6015,
                 "random_R2": 0.7319, "scaffold_R2": 0.6485},
}

for label, s in [("median", sigma_med), ("mean", sigma_mean), ("q75", sigma_q75)]:
    r2_cap = 1 - s * s / var_y
    mae_floor = s * np.sqrt(2 / np.pi)
    out["ceilings"][label] = {
        "sigma": round(s, 4),
        "R2_cap": round(r2_cap, 4),
        "MAE_floor": round(float(mae_floor), 4),
        "random_MAE_over_floor": round(0.5095 / float(mae_floor), 3),
        "scaffold_MAE_over_floor": round(0.6015 / float(mae_floor), 3),
        "random_R2_share_of_cap": round(0.7319 / r2_cap, 3) if r2_cap > 0 else None,
        "scaffold_R2_share_of_cap": round(0.6485 / r2_cap, 3) if r2_cap > 0 else None,
    }

# 另一个视角：多测分子内部，95% 区间有多宽（2 sigma）
out["ci95_width"] = round(2 * 1.96 * sigma_med, 3)

json.dump(out, open(os.path.join(RES, "e11_noise_floor.json"), "w"),
          ensure_ascii=False, indent=2)
print(json.dumps(out, ensure_ascii=False, indent=2))
