#!/usr/bin/env python
"""e9 · 核心对照：同一份数据、同一个模型，只换划分方式

输入
  aidd-mvp/data/processed/dataset_with_scaffold.csv   10468 分子（pIC50 中位数、骨架）
  aidd-mvp/data/processed/X_ecfp.npy                  ECFP4 2048 bit
  aidd-mvp/data/processed/split_indices.json          random / scaffold 两套索引
  aidd_practice/diffdock_demo/results/d2_vina_vs_activity.csv   92 分子（含 Vina 分数）

输出
  results/e9_split_bench.json
  results/e9_predictions.csv        逐分子预测（两种划分的测试集）
  results/e9_domain.csv             测试分子 -> 最近训练集 Tanimoto、绝对误差
  results/e9_92loop.csv             92 分子：真值 / Vina / QSAR / 秩融合

注意
  - 只做回归 pIC50，不涉及分类与 ADMET
  - 预算召回的目标集合 = 当前测试集中真实 pIC50 最高的 30 个（回顾性）
  - 融合分数只用 predict 与 vina_score，绝不引入真值列
"""
import os
import json
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, r2_score
from scipy.stats import spearmanr, kendalltau

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
MVP = os.path.join(ROOT, "aidd-mvp")
PROC = os.path.join(MVP, "data", "processed")
D2 = os.path.join(ROOT, "aidd_practice", "diffdock_demo", "results", "d2_vina_vs_activity.csv")
RES = os.path.join(HERE, "results")
os.makedirs(RES, exist_ok=True)

SEED = 42
N_TREES = 500
TOPFORBUDGET = 30


def log(*a):
    print(*a, flush=True)


def concordance_index(y_true, y_pred):
    """可比较对的排序一致率（忽略并列）。"""
    y_true = np.asarray(y_true, float)
    y_pred = np.asarray(y_pred, float)
    n = len(y_true)
    num = den = 0
    for i in range(n - 1):
        dt = y_true[i] - y_true[i + 1:]
        dp = y_pred[i] - y_pred[i + 1:]
        m = dt != 0
        num += np.sum((dt[m] * dp[m]) > 0)
        den += m.sum()
    return num / den if den else float("nan")


def enrichment_factor(y_true, y_pred, frac=0.05, act_thr=None):
    """EF@x%：前 x% 预测里活性分子的比例 / 全集活性比例。活性阈值取全集 Top-30。"""
    y_true = np.asarray(y_true, float)
    y_pred = np.asarray(y_pred, float)
    n = len(y_true)
    if act_thr is None:
        act_thr = np.sort(y_true)[-TOPFORBUDGET]
    is_act = y_true >= act_thr
    k = max(1, int(round(n * frac)))
    order = np.argsort(-y_pred)
    top = order[:k]
    base = is_act.mean()
    return float(is_act[top].mean() / base) if base > 0 else float("nan"), k


def budget_recall(y_true, y_pred, budgets):
    """目标集 = 真实 Top-30；返回各预算命中数。"""
    y_true = np.asarray(y_true, float)
    y_pred = np.asarray(y_pred, float)
    n = len(y_true)
    target = set(np.argsort(-y_true)[:TOPFORBUDGET])
    order = np.argsort(-y_pred)
    out, expect = {}, {}
    for b in budgets:
        out[b] = int(len(set(order[:b]) & target))
        expect[b] = round(b * TOPFORBUDGET / n, 2)
    return out, expect


def main():
    df = pd.read_csv(os.path.join(PROC, "dataset_with_scaffold.csv"))
    X = np.load(os.path.join(PROC, "X_ecfp.npy"))
    splits = json.load(open(os.path.join(PROC, "split_indices.json")))
    y = df["pIC50_median"].values
    log(f"数据 {X.shape}，pIC50 [{y.min():.2f}, {y.max():.2f}]")

    # 位向量（用于 Tanimoto）
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = []
    for smi in df["smiles"]:
        m = Chem.MolFromSmiles(smi)
        fps.append(gen.GetFingerprint(m) if m is not None else None)

    out = {"n": int(len(df)), "n_features": int(X.shape[1]),
           "n_trees": N_TREES, "seed": SEED, "splits": {}}
    pred_rows, domain_rows = [], []

    for name in ["random", "scaffold"]:
        tr = np.array(splits[name]["train"])
        te = np.array(splits[name]["test"])
        rf = RandomForestRegressor(n_estimators=N_TREES, random_state=SEED, n_jobs=-1)
        rf.fit(X[tr], y[tr])
        p_tr = rf.predict(X[tr])
        p_te = rf.predict(X[te])

        ytr, yte = y[tr], y[te]
        mae = float(mean_absolute_error(yte, p_te))
        rmse = float(np.sqrt(((yte - p_te) ** 2).mean()))
        r2 = float(r2_score(yte, p_te))
        r2_tr = float(r2_score(ytr, p_tr))
        mae_tr = float(mean_absolute_error(ytr, p_tr))
        rho = float(spearmanr(yte, p_te).statistic)
        tau = float(kendalltau(yte, p_te).statistic)
        ci = float(concordance_index(yte, p_te))
        ef5, k5 = enrichment_factor(yte, p_te, 0.05)
        ef10, k10 = enrichment_factor(yte, p_te, 0.10)
        rec, exp = budget_recall(yte, p_te, [10, 20, 30, 100, 200])

        # 骨架重叠（按数据集里的 Murcko 骨架列）
        sc_tr = set(df["scaffold"].values[tr])
        sc_te = df["scaffold"].values[te]
        overlap = int(sum(1 for s in sc_te if s in sc_tr))

        # 最近训练集 Tanimoto
        tr_fps = [fps[i] for i in tr if fps[i] is not None]
        maxsim = []
        for i in te:
            f = fps[i]
            if f is None:
                maxsim.append(0.0)
                continue
            sims = DataStructs.BulkTanimotoSimilarity(f, tr_fps)
            maxsim.append(float(max(sims)))
        maxsim = np.array(maxsim)
        err = np.abs(yte - p_te)

        out["splits"][name] = {
            "n_train": int(len(tr)), "n_test": int(len(te)),
            "train_MAE": round(mae_tr, 4), "train_R2": round(r2_tr, 4),
            "test_MAE": round(mae, 4), "test_RMSE": round(rmse, 4), "test_R2": round(r2, 4),
            "gap_R2": round(r2_tr - r2, 4),
            "spearman": round(rho, 4), "kendall": round(tau, 4), "CI": round(ci, 4),
            "EF5": round(ef5, 3), "EF10": round(ef10, 3), "EF5_k": k5,
            "recall": {str(k): v for k, v in rec.items()},
            "recall_expect": {str(k): v for k, v in exp.items()},
            "scaffold_overlap": overlap,
            "scaffold_overlap_frac": round(overlap / len(te), 4),
            "maxsim": {
                "median": round(float(np.median(maxsim)), 4),
                "q25": round(float(np.quantile(maxsim, .25)), 4),
                "q75": round(float(np.quantile(maxsim, .75)), 4),
                "frac_gt_0.9": round(float((maxsim > 0.9).mean()), 4),
                "frac_gt_0.7": round(float((maxsim > 0.7).mean()), 4),
                "frac_lt_0.4": round(float((maxsim < 0.4).mean()), 4),
            },
            "MAE_median_abs": round(float(np.median(err)), 4),
            "MAE_p90": round(float(np.quantile(err, .9)), 4),
        }
        log(f"[{name}] train R2={r2_tr:.4f}  test MAE={mae:.4f} R2={r2:.4f} "
            f"rho={rho:.4f} CI={ci:.4f} EF5={ef5:.2f} recall={rec} maxsim_med={np.median(maxsim):.3f}")

        for j, i in enumerate(te):
            pred_rows.append({"split": name, "smiles": df["smiles"].values[i],
                              "scaffold": df["scaffold"].values[i],
                              "true_pIC50": float(yte[j]), "pred_pIC50": float(p_te[j]),
                              "abs_err": float(err[j]), "maxsim": float(maxsim[j])})
        for j in range(len(te)):
            domain_rows.append({"split": name, "maxsim": float(maxsim[j]),
                                "abs_err": float(err[j]), "true": float(yte[j]),
                                "pred": float(p_te[j])})

        # 适用范围：域内 / 域外（阈值 0.4，与后文口径一致）
        for thr in (0.7, 0.5, 0.4):
            inm = maxsim >= thr
            if inm.sum() > 5 and (~inm).sum() > 5:
                out["splits"][name][f"domain_R2_{thr}"] = {
                    "in_n": int(inm.sum()), "out_n": int((~inm).sum()),
                    "in_R2": round(float(r2_score(yte[inm], p_te[inm])), 4),
                    "out_R2": round(float(r2_score(yte[~inm], p_te[~inm])), 4),
                    "in_MAE": round(float(np.abs(yte[inm] - p_te[inm]).mean()), 4),
                    "out_MAE": round(float(np.abs(yte[~inm] - p_te[~inm]).mean()), 4),
                }

    # ---------- 92 分子闭环：Vina / QSAR / 秩融合 ----------
    d2 = pd.read_csv(D2)
    log(f"92 分子表 {d2.shape}")
    keep = ~df["smiles"].isin(set(d2["smiles"]))
    rf92 = RandomForestRegressor(n_estimators=N_TREES, random_state=SEED, n_jobs=-1)
    rf92.fit(X[keep.values], y[keep.values])
    idx_map = {s: i for i, s in enumerate(df["smiles"].values)}
    qsar = np.array([rf92.predict(X[[idx_map[s]]])[0] for s in d2["smiles"]])
    yt = d2["pIC50_median"].values
    vina = d2["vina_score"].values

    # 秩融合：Vina 越低越好 -> 取负；QSAR 越高越好。秩平均，不用任何真值。
    r_v = (-vina).argsort().argsort()
    r_q = qsar.argsort().argsort()
    fused = 0.5 * (r_v + r_q)

    rows = []
    for label, score in [("Vina", -vina), ("QSAR", qsar), ("RankFusion", fused)]:
        rec92, exp92 = budget_recall(yt, score, [10, 20, 30])
        rows.append({"strategy": label,
                     "recall_10": rec92[10], "recall_20": rec92[20], "recall_30": rec92[30],
                     "spearman": round(float(spearmanr(yt, score).statistic), 4),
                     "MAE": round(float(np.abs(yt - score).mean()), 4) if label != "Vina" else None})
        log(f"  92loop {label:10s} recall={rec92} rho={rows[-1]['spearman']}")
    out["loop92"] = {"n": int(len(d2)), "rows": rows,
                     "expect": {str(k): v for k, v in exp92.items()},
                     "note": "融合只用 Vina 与 QSAR 预测，未引入 pIC50 真值列"}

    pd.DataFrame(pred_rows).to_csv(os.path.join(RES, "e9_predictions.csv"), index=False)
    pd.DataFrame(domain_rows).to_csv(os.path.join(RES, "e9_domain.csv"), index=False)
    pd.DataFrame({"smiles": d2["smiles"], "true": yt, "vina": vina, "qsar": qsar,
                  "fused": fused}).to_csv(os.path.join(RES, "e9_92loop.csv"), index=False)
    json.dump(out, open(os.path.join(RES, "e9_split_bench.json"), "w"),
              ensure_ascii=False, indent=2)
    log("wrote results/e9_split_bench.json")


if __name__ == "__main__":
    main()
